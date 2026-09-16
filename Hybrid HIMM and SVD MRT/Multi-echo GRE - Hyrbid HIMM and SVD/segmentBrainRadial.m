function [skullMask, bodyMask, brainMask] = segmentBrainRadial(I, body_mask_mGRE)
% segmentBrainRadial: Create skull, body, and brain masks for HIMMmotion.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to the license file.
%
% Inputs:
%   I               : Gradient-echo image of a transverse brain slice
%   body_mask_mGRE  : Initial tissue/body mask
%
% Outputs:
%   skullMask       : Mask of the outer skull region
%   bodyMask        : Mask enclosed by the detected inner skull boundary
%   brainMask       : Brain mask obtained by eroding bodyMask
%
% Requires Image Processing Toolbox.

    %% Normalize and smooth image

    I = double(abs(I(:,:,1)));
    I = I - min(I(:));

    if max(I(:)) > 0
        I = I / max(I(:));
    end

    I_smooth = imgaussfilt(I, 1.0);


    %% Determine outer head mask

    threshold = graythresh(I_smooth);

    headMask = I_smooth > 0.7 * threshold;

    headMask = imclose(headMask, strel('disk', 3));
    headMask = imfill(headMask, 'holes');
    headMask = bwareaopen(headMask, 200);
    headMask = largestComponent(headMask);

    % Additional closing to bridge gaps in the outer head boundary
    headMask = imclose(headMask, strel('disk', 5));
    headMask = imfill(headMask, 'holes');
    headMask = largestComponent(headMask);


    %% Determine center and approximate radius of the head

    stats = regionprops(headMask, 'Centroid', 'EquivDiameter');

    cx = stats.Centroid(1);
    cy = stats.Centroid(2);

    headRadius = stats.EquivDiameter / 2;


    %% Radial search parameters

    nAngles = 360;

    theta = linspace(0, 2*pi, nAngles + 1);
    theta(end) = [];

    maxR = ceil(1.2 * headRadius);

    % Search region inward from the outer head boundary
    inwardMin = max(4,  round(0.10 * headRadius));
    inwardMax = max(10, round(0.25 * headRadius));

    innerRadius = nan(1, nAngles);
    outerRadius = nan(1, nAngles);


    %% Detect dark gap between brain and skull

    for k = 1:nAngles

        th = theta(k);
        r = 0:maxR;

        x = cx + r .* cos(th);
        y = cy + r .* sin(th);

        % Sample intensity and head mask along radial line
        intensityProfile = interp2(I_smooth, x, y, 'linear', 0);
        maskProfile = interp2(double(headMask), x, y, 'nearest', 0) > 0.5;

        % Find outer head boundary
        idxOuter = find(maskProfile, 1, 'last');

        if isempty(idxOuter) || idxOuter < inwardMin + 3
            continue;
        end

        outerRadius(k) = r(idxOuter);

        % Search for dark gap inside outer head boundary
        idxStart = max(2, idxOuter - inwardMax);
        idxEnd   = idxOuter - inwardMin;

        if idxEnd <= idxStart
            continue;
        end

        searchProfile = intensityProfile(idxStart:idxEnd);
        searchProfile = smoothdata(searchProfile, 'gaussian', 7);

        % Find darkest point
        [~, idxMin] = min(searchProfile);
        idxGap = idxStart + idxMin - 1;

        % Refine using nearby local minimum
        localIdx = localMinimumNear( ...
            intensityProfile, ...
            idxGap, ...
            max(idxStart, idxGap - 3), ...
            min(idxEnd, idxGap + 3));

        if ~isempty(localIdx)
            idxGap = localIdx;
        end

        innerRadius(k) = r(idxGap);
    end


    %% Fill missing values and smooth detected boundary

    innerRadius = fillmissingCircular(innerRadius);
    outerRadius = fillmissingCircular(outerRadius);

    % Ensure detected inner boundary remains within the head
    innerRadius = min(innerRadius, outerRadius);
    innerRadius(innerRadius < 1) = 1;

    % Smooth contour across neighboring angles
    innerRadius = smoothCircular(innerRadius, 21);
    innerRadius = innerRadius - 5;


    %% Convert radial contour to body mask

    xBoundary = cx + innerRadius .* cos(theta);
    yBoundary = cy + innerRadius .* sin(theta);

    bodyMask = poly2mask( ...
        xBoundary, ...
        yBoundary, ...
        size(I,1), ...
        size(I,2));

    bodyMask = bodyMask & headMask;

    % Clean mask
    bodyMask = imfill(bodyMask, 'holes');
    bodyMask = imopen(bodyMask, strel('disk', 2));
    bodyMask = imclose(bodyMask, strel('disk', 4));
    bodyMask = largestComponent(bodyMask);


    %% Create brain and skull masks

    erosionRadius = 1;
    brainMask = imerode(bodyMask, strel('disk', erosionRadius));

    % Region between bodyMask and brainMask
    skullMask = bodyMask & ~brainMask;

    % Only retain voxels included in the original mGRE body mask
    skullMask = skullMask & logical(body_mask_mGRE);

end


function BW = largestComponent(BW)
% Keep only the largest connected component.

    CC = bwconncomp(BW);

    if CC.NumObjects == 0
        BW = false(size(BW));
        return;
    end

    componentSize = cellfun(@numel, CC.PixelIdxList);
    [~, idxLargest] = max(componentSize);

    BW = false(size(BW));
    BW(CC.PixelIdxList{idxLargest}) = true;

end


function idx = localMinimumNear(profile, idx0, idxStart, idxEnd)
% Find a local minimum close to the initially detected minimum.

    idx = [];

    for k = max(idxStart + 1, idx0 - 2):min(idxEnd - 1, idx0 + 2)

        if profile(k) <= profile(k - 1) && ...
           profile(k) <= profile(k + 1)

            idx = k;
            return;

        end
    end

end


function y = smoothCircular(x, windowSize)
% Circular moving-average smoothing.

    if mod(windowSize, 2) == 0
        windowSize = windowSize + 1;
    end

    halfWindow = floor(windowSize / 2);

    xPadded = [ ...
        x(end-halfWindow+1:end), ...
        x, ...
        x(1:halfWindow)];

    kernel = ones(1, windowSize) / windowSize;

    yPadded = conv(xPadded, kernel, 'same');

    y = yPadded(halfWindow+1 : halfWindow+numel(x));

end


function x = fillmissingCircular(x)
% Fill missing values in a circular 1D signal using linear interpolation.

    n = numel(x);
    valid = ~isnan(x);

    if all(valid)
        return;
    end

    if ~any(valid)
        x(:) = 1;
        return;
    end

    idx = 1:n;

    xExtended = [x, x, x];
    idxExtended = [idx, idx+n, idx+2*n];

    validExtended = ~isnan(xExtended);

    x = interp1( ...
        idxExtended(validExtended), ...
        xExtended(validExtended), ...
        idx+n, ...
        'linear');

end