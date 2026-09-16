function [m,wts,Ac,c,A,phi,phase_motion,c_motion] = OptimizeHIMMSVD(imgs, algp, scanp, lib, Acswitch,Acmotionswitch,mswitch,phiswitch,minit,phiinit)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% OptimizeHIMMSVD: main function call
%
% Function to solve for heating, field shifts, select a baseline, and 
% correct for motion. Note that the algorithm only works if the region of
% interest is surrounded by a layer of fat
%
% Algorithm code for the iterative solver was dfeveloped by Megan Poorman 
% and William Grissom. The code presented in this version does not use the 
% regularization term anymore as it is not focused on applications using
% localized heating. Alternatively, the drift field (Ac) was initialized 
% using near-harmonic 2D reconstruction (proposed by Salomir et al.).
% 
% Furthermore, the algorithm is also extended by a motion-correction
% scheme. Head motion was found to affect the phase as well. Luckily, the
% change in the B0 field due to motion is linear with respect to the amount
% of rotation/translation. If we know how the head has moved during dynamic
% scanning and we know how the B0 field behaves with respect to certain
% rotations/translations, we can compensate for motion-induced B0
% fluctuations.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x Necho)
%   ----- algorithm params
%     algp.B0_motion = subject-specific motion fields;
%     algp.order = order of the spatial function used to fit drift;
%     algp.method = functionality of approaching drift fields (polymonials or spherical harmonics);
%     algp.nmiter = number of temperature (m) fit iterations;
%     algp.nciter = number of drift (c) fit iterations;
%     algp.Acswitch = whether drift fields will be fitted;
%     algp.Acmotionswitch = whether motion fields will be fitted;
%     algp.positivity_constrain = constrain temperature to be positive whilst trying to solve for drift and motion fields;
%     algp.max_iters = maximum number of iterations;
%     algp.stopthresh = fit tolerance;
%     algp.suppress_info = whether updates will be visualized on the command window;
%     algp.rad2degC = conversion factor between radians and degree celsius;
%
%   ----- masking parameters
%     algp.body_factor = value used to vary the threshold for the in-phase image;
%     algp.Nstd_ol = NH2DRecon: filter out local extreme values in the fat border;
%     algp.neighbours = NH2DRecon: define local neighbourhood;
%     algp.min_nb = NH2DRecon: eliminate voxels at the fat-tissue interface from the fat mask;
%
%   ----- scan parameters
%     scanp.dim = [rows cols];
%     scanp.b0 = b0;
%     scanp.tes = tes;
%     scanp.prc = prc;
%     scanp.T_rigid = motion parameters for the dynamic being evaluated;
% 
%   ----- baseline library
%     lib.Wlib = Wlib;
%     lib.Flib = Flib;
%     lib.dw0lib = dw0lib;
%     lib.R2starlib = R2starlib;
%     lib.fatmodel = fatmodel;
%
%   Acswitch:           integer specifying whether we should include a fit
%                       for drift in the algorithm
%   Acmotionswitch:     integer specifying whether we should include a fit
%                       for motion in the algorithm
%   mswitch:            boolean whether we should include a fit
%                       for temperature in the algorithm
%   phiswitch:          boolean whether we should include a fit
%                       for Tx/Rx gain in the algorithm
%   minit:              initialization for temperature
%   phiinit:            initialization for Tx/Rx gain
%
%
% Output: 
%   m:                  frequency shift due to heating
%   wts:                baseline weight
%   Ac:                 estimated drift field
%   c:                  spherical harmonics coefficients fit
%   A:                  spherical harmonics basis functions
%   phi:                Tx/Rx gain estimated
%   phase_motion:       estimated motion field
%   c_motion:           motion coefficients fit


%% --- Initialization --- %%

%--- Define scan size
if isscalar(scanp.dim) %assume square
    scanp.dim = [scanp.dim scanp.dim];
end

%--- Initialize universal DC offsets
if ~exist('phiinit','var')
    phi = zeros(scanp.dim(1),scanp.dim(2));
else
    phi = phiinit;
end 

%--- Initialize field drifts
if strcmp(algp.method,'poly')
    A = get_polynomial_function(scanp.dim,algp.order);
elseif strcmp(algp.method,'spher')
    A = get_spherical_function(scanp.dim,algp.order);
else
    message = {'The function with coefficients that are to be '
        'optimized was specified incorrectly, '
        'please set algp.method to the correct string. '};
    errorStruct.message = sprintf('%s\n',message{:});
    errorStruct.identifier = 'MyFunction:parameterIncorrecctlyDefined';
    error(errorStruct)
end
c = zeros(size(A,2),1);
if scanp.dim(1)==scanp.dim(2)
    Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
else
    Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
end
Ac_salomir = Ac;


%--- Preparation for motion fitting
% Initialize the motion matrix A that includes changes in B0 rotations and translations
A_motion = [];
for motion_ix = 1:size(algp.B0_motion,3)
    motion_function = algp.B0_motion(:,:,motion_ix);
    A_motion = [ A_motion motion_function(:) ];
end
c_motion = zeros(size(A_motion,2),1);
phase_motion = zeros(size(Ac));

%--- Get masks
[fat_mask,body_mask,brain_mask] = segmentBrainRadial(scanp.dyn_img_mask, algp.body_mask);



%% --- Start HIMM iterations with motion compensation by means of SVD --- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Summary of the current algorithm code

% The code consists of 4 separate sections
% 1. Finding the best-matching reference
%    The library of references contains water/fat reconstructions of several
%    dynamics acquired in the 'baseline' scenario. A constrained least-squares
%    problem is solved that minimizes the difference between the measured 
%    dynamic and a weighted combination of reference dynamics (updated for 
%    a change in drift field), subject to the constraint that the coefficients 
%    sum to one.
%
% 2. Near Harmonic 2D reconstruction (initialization phase)
%    Motion correction is performed based on the measured motion
%    parameters. Subsequently, the motion-corrected phase accumulation
%    (between the current dynamic and the determined reference dynamic) is
%    determined. This phase accumulation is used for fat-referenced
%    near-harmomic 2D reconstruction.
%
% 3. Optimization (with positivity constrain), optional
%    A jointly optimization of motion fields, drift fields, and temperature 
%    fields is performed. Here, temperature changes are constrained to be 
%    positive such that the motion and drift fields converge in the right 
%    direction. 
%
% 4. Optimization (without positivity constrain)
%    A jointly optimization of motion fields, drift fields, and temperature 
%    fields is performed. Here, temperature changes are not constrained to 
%    be positive.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization factor (median of heating dynamics * num echoes)
algp.norm = abs(median(imgs,'all') * length(scanp.tes));

% Create B0 mask (region in which motion correction can be applied)
mask_B0 = body_mask;
mask_B0(sum(algp.B0_motion,3) == 0) = 0;

% Set masks for fitting
fit_mask_drift = brain_mask;   % drift field fitting mask
fit_mask_motion = brain_mask;  % motion field fitting mask

m = minit.*brain_mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Find the combination of dynamics in the reference library that 
%    best represents the current dynamic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if algp.suppress_info ~= 1
    fprintf('    Finding best-matching reference\n');
end

% Build Z matrix using previously estimated values of Ac to correct for the
% change in drift field
for jj = 1:size(lib.Wlib,3)

    % Individually construct every model-based reference dynamic 
    % The drift field will be updated in each reference dynamic, and each
    % reference will have a slightly different motion field
    wts = zeros([size(lib.Wlib,3) 1])/size(lib.Wlib,3);
    wts(jj) = 1;

    if Acmotionswitch
        % Initialize c_motion with rotation parameters corresponding to the
        % measured values of the registration (scanp.T_rigid)
        c_motion = zeros(size(A_motion,2),1);
        for i_type = 1:size(A_motion,2)
            for i_ref = 1:size(scanp.T_rigid,2)
                c_motion(i_type) = c_motion(i_type) + wts(i_ref)*scanp.T_rigid(i_type,i_ref);
            end
        end
    
        % Calculate motion induced field inhomogeneity map
        if scanp.dim(1)==scanp.dim(2)
            phase_motion = reshape(A_motion*c_motion,[scanp.dim(1) scanp.dim(2)]);
        else
            phase_motion = reshape(A_motion*c_motion,[scanp.dim(2) scanp.dim(1)]).';
        end
    end

    % Compute model values
    Z(:,:,:,jj) = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                             scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ... 
                             lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);

end
Z = Z.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 
wts = f_update_lsqlin(imgs,Z,brain_mask);

% Update normalization factor: median of heating dynamics * num echoes
body_masks = repmat(body_mask,1,1,length(scanp.tes),size(imgs,4));
imgsdyn_masked = imgs; imgsdyn_masked(body_masks ~= 1) = nan; 
imgsdyn_masked = imgsdyn_masked(:);
imgsdyn_masked = imgsdyn_masked(~isnan(imgsdyn_masked));
algp.norm = abs(median(imgsdyn_masked) * length(scanp.tes));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. initialize and perform near-harmonic 2D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize motion field with current m and Ac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Acmotionswitch
    
    % Initialize c_motion with rotation parameters corresponding to the
    % measured values of the registration (scanp.T_rigid)
    c_motion = zeros(size(A_motion,2),1);
    for i_type = 1:size(A_motion,2)
        for i_ref = 1:size(scanp.T_rigid,2)
            c_motion(i_type) = c_motion(i_type) + wts(i_ref)*scanp.T_rigid(i_type,i_ref);
        end
    end

    % Calculate motion induced field inhomogeneity map
    if scanp.dim(1)==scanp.dim(2)
        phase_motion = reshape(A_motion*c_motion,[scanp.dim(1) scanp.dim(2)]);
    else
        phase_motion = reshape(A_motion*c_motion,[scanp.dim(2) scanp.dim(1)]).';
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform near-harmonic 2D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Acswitch == 1

    fit_mask = fat_mask;
    
    for kk = 1:algp.nmiter

        % Calculate and apply weights to baseline images
        Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);

        % We want to get a first estimate of Ac by performing near-harmonic 2D
        % reconstruction on the fat. For that, we first need to calculate the
        % phase accumulation in fat
        for jj = 1:size(lib.Wlib,3)
            % Compute model values
            Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
             scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
             lib.R2starlib(:,:,jj),0,scanp.prc,lib.fatmodel);
            Ztot = Ztot + Z*wts(jj);
        end
        Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 

        % Update drift field using the fit mask (fat-only)
        c = c_update(imgs(:,:,1:end),Ztot(:,:,1:end),A,c,scanp.tes(1:end),scanp.prc,fit_mask);

        % Recalculate frequency shift map
        if scanp.dim(1)==scanp.dim(2)
            Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
        else
            Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
        end
        
    end

    %%%%%%%%%%%%%
    % run Salomir
    %%%%%%%%%%%%%
    
    % fprintf('    Performing Near-Harmonic 2D reconstruction\n');
    [Ac,c,Ac_salomir] = NH2DReconExecution(A,Ac,Ac_salomir,fit_mask,body_mask,brain_mask,algp,scanp);

else

    Ac = zeros(size(Ac));
    Ac_salomir = zeros(size(Ac));
    c = zeros(size(c));
    
end

% Set initial temperature
% Determine using near-harmonic 2D reconstruction (and a positivity constraint for initialization-only)
if mswitch
    for kk = 1:algp.nmiter

        % Build Z matrix using true values of Ac and m
        Zw = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
        Zf = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
        
        for jj = 1:size(lib.Wlib,3)
            % Compute model values of water and fat separately
            [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ... 
                lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
            Zw = Zw + foow*wts(jj);
            Zf = Zf + foof*wts(jj);
        end

        % Update m
        m = m_update(imgs(:,:,1:end)-Zf(:,:,1:end),Zw(:,:,1:end),m,scanp.tes(1:end),scanp.prc,mask_B0);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run phi_update with current f, m, Ac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This accounts for any universal DC offsets due to reciever gain etc.
if phiswitch
    for kk = 1:algp.nciter
        
        % Calculate and apply weights to baseline images
        Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);

        for jj = 1:size(lib.Wlib,3)
            % Compute model values
            Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
            Ztot = Ztot + Z*wts(jj);
        end
        Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 

        % Update phi
        phi = phi_update(imgs,Ztot,phi(:),scanp.tes);

        % Recalculate frequency shift map
        phi = reshape(phi,[scanp.dim(1) scanp.dim(2)]);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. start the iterative process of fitting with spherical harmonics
%    Executed with a positivity constrain in temperature
%    Ac is initialized properly through near-harmonic 2D reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if algp.positivity_constrain == 1

    % calculate a terrible cost that is guaranteed to get iterations going
    cost = cost_eval(imgs,scanp,lib,lib.dw0lib,phase_motion,Ac,m,phi,algp,ones(size(imgs(:,:,1)))); 
    cost = 2*cost;
    ii = 0; % iteration counter
    
    % The optimization algorithm is run for at least 1 iterations
    converged = false;
    while converged == false || ii <= 0
    
        iter_stop = algp.max_iters - 1; % maximum of iterations
        if ii > iter_stop
            break
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run c_update with current f, m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Acswitch
            for kk = 1:algp.nciter
            
                % Calculate and apply weights to baseline images
                Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
                
                for jj = 1:size(lib.Wlib,3)
                    % Compute model values
                    Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                        lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                    Ztot = Ztot + Z*wts(jj);
                end
                Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3)]);
                
                % Update c
                c = c_update(imgs(:,:,1:end),Ztot(:,:,1:end),A,c,scanp.tes(1:end),scanp.prc,fit_mask_drift);
                
                % Recalculate frequency shift map
                if scanp.dim(1)==scanp.dim(2)
                    Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
                else
                    Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run c_motion_update with current m and Ac
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if Acmotionswitch
    
            % Use gradient descent to find coefficients for motion fields
            for kk = 1:algp.nciter
    
                % Use gradient descent to find correct weights for motion fields
                Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
                
                for jj = 1:size(lib.Wlib,3)
                    % Compute model values
                    Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                        lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                    Ztot = Ztot + Z*wts(jj);
                end
                Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3)]);
    
                % Update c
                lambda = 1.0;
                c_motion = c_motion_update(imgs(:,:,1:end),Ztot(:,:,1:end),A_motion,c_motion,scanp.tes(1:end),scanp.prc,fit_mask_motion,lambda);
                
                if scanp.dim(1)==scanp.dim(2)
                    phase_motion = reshape(A_motion*c_motion,[scanp.dim(1) scanp.dim(2)]);
                else
                    phase_motion = reshape(A_motion*c_motion,[scanp.dim(2) scanp.dim(1)]).';
                end
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run phi_update with current f, m, Ac
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This accounts for any universal DC offsets due to reciever gain etc.
        if phiswitch
            for kk = 1:algp.nciter
                
                % Calculate and apply weights to baseline images
                Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
     
                for jj = 1:size(lib.Wlib,3)
                    % Compute model values
                    Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                        lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                    Ztot = Ztot + Z*wts(jj);
                end
                Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 
    
                % Update phi
                phi = phi_update(imgs,Ztot,phi(:),scanp.tes);
    
                % Recalculate frequency shift map
                phi = reshape(phi,[scanp.dim(1) scanp.dim(2)]);
    
            end
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run m_update with current f, Ac
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if mswitch
            for kk = 1:algp.nmiter
    
                % Build Z matrix using true values of Ac and m
                Zw = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
                Zf = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
                
                for jj = 1:size(lib.Wlib,3)
                    % Compute model values of water and fat separately
                    [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                        scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ... 
                        lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                    Zw = Zw + foow*wts(jj);
                    Zf = Zf + foof*wts(jj);
                end
    
                % Update m
                m = m_update(imgs(:,:,1:end)-Zf(:,:,1:end),Zw(:,:,1:end),m,scanp.tes(1:end),scanp.prc,mask_B0);
                m(brain_mask == 0) = 0;
                m((m > 0) & (mask_B0 == 1)) = 0;

            end
        end
        
        % report cost + iteration
        costOld = cost; ii = ii + 1;
        cost = cost_eval(imgs,scanp,lib,wts,phase_motion,Ac,m,phi,algp,brain_mask);
        converged = costOld - cost < algp.stopthresh*costOld;
    
        if algp.suppress_info ~= 1
            fprintf('    Optimization iteration %d (positivity constrain in dT): Cost = %0.2d | difference = %0.2d (%d)\n',ii,cost,cost-costOld,converged);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. start the iterative process without positivity constrain in dT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate a terrible cost that is guaranteed to get iterations going
cost = cost_eval(imgs,scanp,lib,wts,phase_motion,Ac,m,phi,algp,ones(size(imgs(:,:,1)))); 
cost = 2*cost;
ii = 0; % iteration counter

% The optimization algorithm is run for at least 1 iteration
converged = false;
while converged == false || ii <= 0

    iter_stop = algp.max_iters - 1; % maximum of iterations
    if ii > iter_stop
        break
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run c_update with current f, m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Acswitch
        for kk = 1:algp.nciter
        
            % Calculate and apply weights to baseline images
            Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            
            for jj = 1:size(lib.Wlib,3)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3)]);
                        
            % Update c
            c = c_update(imgs(:,:,1:end),Ztot(:,:,1:end),A,c,scanp.tes(1:end),scanp.prc,fit_mask_drift);
            
            % Recalculate frequency shift map
            if scanp.dim(1)==scanp.dim(2)
                Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
            else
                Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
            end
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run c_motion_update with current m and Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Acmotionswitch

        % Use gradient descent to find coefficients for motion fields
        for kk = 1:algp.nciter

            % Use gradient descent to find correct weights for motion fields
            Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            
            
            for jj = 1:size(lib.Wlib,3)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3)]);
            
            % Update c
            lambda = 1.0;
            c_motion = c_motion_update(imgs(:,:,1:end),Ztot(:,:,1:end),A_motion,c_motion,scanp.tes(1:end),scanp.prc,fit_mask_motion,lambda);
            
            if scanp.dim(1)==scanp.dim(2)
                phase_motion = reshape(A_motion*c_motion,[scanp.dim(1) scanp.dim(2)]);
            else
                phase_motion = reshape(A_motion*c_motion,[scanp.dim(2) scanp.dim(1)]).';
            end
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run phi_update with current f, m, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This accounts for any universal DC offsets due to reciever gain etc.
    if phiswitch
        for kk = 1:algp.nciter
            
            % Calculate and apply weights to baseline images
            Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
 
            for jj = 1:size(lib.Wlib,3)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 

            % Update phi
            phi = phi_update(imgs,Ztot,phi(:),scanp.tes);

            % Recalculate frequency shift map
            phi = reshape(phi,[scanp.dim(1) scanp.dim(2)]);

        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run m_update with current f, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mswitch
        for kk = 1:algp.nmiter

            % Build Z matrix using true values of Ac and m
            Zw = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            Zf = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            
            for jj = 1:size(lib.Wlib,3)
                % Compute model values of water and fat separately
                [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ... 
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Zw = Zw + foow*wts(jj);
                Zf = Zf + foof*wts(jj);
            end

            % Update m
            m = m_update(imgs(:,:,1:end)-Zf(:,:,1:end),Zw(:,:,1:end),m,scanp.tes(1:end),scanp.prc,mask_B0);

        end
    end
    
    % report cost + iteration
    costOld = cost;ii = ii + 1;
    cost = cost_eval(imgs,scanp,lib,wts,phase_motion,Ac,m,phi,algp,brain_mask);
    converged = costOld - cost < algp.stopthresh*costOld;

    if algp.suppress_info ~= 1
        fprintf('    Optimization iteration %d (no positivity constrain in dT): Cost = %0.2d | difference = %0.2d (%d)\n',ii,cost,cost-costOld,converged);
    end
end

% Perform temperature determination again, but now without the masks and negativity constrain
if mswitch
    for kk = 1:algp.nmiter

        % Build Z matrix using true values of Ac and m
        Zw = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
        Zf = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
        
        for jj = 1:size(lib.Wlib,3)
            % Compute model values of water and fat separately
            [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ... 
                lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
            Zw = Zw + foow*wts(jj);
            Zf = Zf + foof*wts(jj);
        end

        % Update m
        m = m_update(imgs(:,:,1:end)-Zf(:,:,1:end),Zw(:,:,1:end),m,scanp.tes(1:end),scanp.prc,ones(size(brain_mask)));

    end
end


%% Cost function
% Function to evaluate the overall cost function
function cost = cost_eval(imgs,scanp,lib,wts,phase_motion,Ac,m,phi,algp,mask)
    Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
    for jj = 1:size(lib.Wlib,3)
        Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
            scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac+phase_motion, ...
            lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
        Ztot = Ztot + Z*wts(jj);
    end

    Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]);
    imgs_input = imgs.*repmat(mask,1,1,length(scanp.tes));
    Ztot_input = Ztot.*repmat(mask,1,1,length(scanp.tes));
    cost = norm((imgs_input(:)-Ztot_input(:))/algp.norm)^2;