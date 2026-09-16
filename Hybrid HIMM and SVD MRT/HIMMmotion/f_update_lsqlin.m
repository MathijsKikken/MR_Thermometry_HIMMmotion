function [w,libinds] = f_update_lsqlin(imgs,Z,mask)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% f_update: find the weights that best match the baseline library to the
% multi-echo images using constrained least squares (lsqlin)
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs:   Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:      Model image for each echo and each baseline
%           (Nx x Ny x Necho x Nbaseline)
%   mask:   Spatial mask (Nx x Ny)
%
% Outputs:
%   w:        baseline weights
%   libinds:  library indices used (here always 1:Nbaseline)

    nBaseline = size(Z,4);

    % Only one baseline image
    if nBaseline == 1
        w = 1;
        libinds = 1;
        return
    end

    % Build 3D mask over all echoes
    mask3 = repmat(mask,[1 1 size(imgs,3)]) ~= 0;

    % Vectorize measured image inside mask
    y = imgs(mask3);
    y = y(:);

    % Build library matrix A: one column per artificial image
    nVox = numel(y);
    A = zeros(nVox, nBaseline, 'like', y);

    for k = 1:nBaseline
        Zk = Z(:,:,:,k);
        A(:,k) = Zk(mask3);
    end

    % Handle complex-valued data by stacking real and imaginary parts
    % lsqlin works with real variables, so rewrite:
    % ||A*w - y||^2  ->  ||[real(A); imag(A)]*w - [real(y); imag(y)]||^2
    if ~isreal(A) || ~isreal(y)
        A = [real(A); imag(A)];
        y = [real(y); imag(y)];
    end

    % Equality constraint: coefficients must sum to 1
    Aeq = ones(1,nBaseline);
    beq = 1;

    % Bounds: convex combination
    lb = zeros(nBaseline,1);
    ub = ones(nBaseline,1);

    % Solve constrained least-squares problem
    opts = optimoptions('lsqlin', ...
        'Algorithm','interior-point', ...
        'Display','off', ...
        'OptimalityTolerance',1e-10, ...
        'StepTolerance',1e-12, ...
        'MaxIterations',1e5);

    [w,~,~,~,~] = lsqlin(A,y,[],[],Aeq,beq,lb,ub,[],opts);

    % Library indices used
    libinds = 1:nBaseline;
end