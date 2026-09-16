%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% 
%  HIMM: Harmonic Initialized Model-based Multi-echo  
%  SVD: Singular Value Decompotion (used as a motion correction scheme)
%  Main script to extract temperature, drift, and motion fields from measured data
%
%  Besides B0 field drift and temperature, the phase over time will also 
%  change due to cardiac, respiratory, and head motion. 
%   - Cardiac motion is accounted for through cardiac triggered. 
%   - Respiratory motion is accounted for through a multi-referenced approach
%   - Head motion is accounted for through a motion-correction scheme
%     This scheme used characteristic-motion induced field maps (every
%     type of motion will result in very characteristic changed in the B0
%     field). This means that head motion-induced field changes can be
%     compensated for if we know how the head moved during dynamic scanning
%     and if we know what these charachteristic motion-induced fiels look
%     like. The latter can be obtained from the following GitHub
%     repository: https://github.com/MathijsKikken/GetMotionMaps
%
%  First, the algorithm determines which of the reference dynamics in the
%  library is most similar to the dynamic image being evaluated.
%  
%  Subsequenly, an initial estimate of the B0 field drift is acquired
%  through near-harmonic 2D reconstruction. The skull is insensitive to
%  temperature-induced phase changes and can therefore be used to provide a
%  good initial estimate of the B0 field drift in the tranverse plane
%  (because the skull completely surrounds the brain in this plane).
%  
%  After this initial estimate of the B0 field drift was acquired. The
%  algorithm jointly updates the motion-induced field disturbances (motion 
%  fields), field drift, and the temperature map to make sure the model 
%  best resembles the dynamic currently being evaluated.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% The basis of this algorithm was created by Megan Poorman & William Grissom
% (https://github.com/poormanme/waterFatSeparated_MRThermometry)
%
% Updated from: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Do not reproduce, distribute, or modify without proper citation according
% to license file

clear all; close all; clc;

% Add the software (HIMM algorithm code and ISMRM toolbox) to the directory
addpath(genpath('.\Code\MR_Thermometry_HIMMmotion\Hybrid HIMM and SVD MRT\HIMMmotion'))
addpath(genpath('.\Code\ISMRM Water Fat Toolbox')) % ISMRM fat/water toolbox


%% Define scanner and sequence parameters

b0 = 7;                     % Tesla
frequency_system = 298*1e6; % Frequency measured by used MR system [Hz]
prc = 1;                    % Direction of precession (scanner-dependent) 
                            % +1 = heating induces negative apparent freq shift
                            % -1 = heating induces positive apparent freq shift
phi0shift = 0;              % Tx/Rx gain
alpha = -0.01;              % ppm/deg C
gyro = 42.5778;             % MHz/T
rad2degC = 1/(2*pi*b0*alpha*gyro);  % Factor to convert radians to temperature


%% Provide subject data and define the scan that will be evaluated

% Subject definition
datafile = '.\Data\';
Pt       = 1;                   % subject number

% Scan definition
sequence         = 'heating';  % 'baseline' or 'heating'
info_name_append = '';          % '' for no appended info, '_motion_sag' for motion scan
slice            = 1;           % usually set to 1, except for M2D acquisitions
num_stacks       = 1;           % Number of transverse images acquired in a single dynamic (usually 1)
discard_dyns     = 10;          % Number of dynamics that are discarded
num_avgs         = 10;           % Number of dynamics in the reference library


%% Load corresponding rigid transformation parameters and determine references

%T_rigid_temperature = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_motion_cor.nii.gz']);
T_rigid_temperature = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_',sequence,info_name_append,'.nii.gz']);

% Provide a visualization of the rotation/translation parameters over time
figure; 
subplot(4,1,1); plot(1:size(T_rigid_temperature,2),T_rigid_temperature(2,:),'-','Color',"#0072BD",LineWidth=3); title('Rotation around the LR axis'); hold on;
subplot(4,1,2); plot(1:size(T_rigid_temperature,2),T_rigid_temperature(3,:),'-','Color',"#FF0000",LineWidth=3); title('Rotation around the AP axis'); hold on;
subplot(4,1,3); plot(1:size(T_rigid_temperature,2),T_rigid_temperature(6,:),'-','Color',"#000000",LineWidth=2); title('FH translations'); hold on;
subplot(4,1,4); plot(1:size(T_rigid_temperature,2),T_rigid_temperature(1,:),'-','Color',"#00FF00",LineWidth=2); title('Rotation around the FH axis'); hold on;


%% Load scan and define which (averaged) dynamics should be used as reference 

% Specify which echoes are used for the reconstruction
relevant_echoes = [1,2:16];
img_name_append = strcat(info_name_append,'_registered_slice',num2str(slice-1));  % '' for no appended info, '_shim2' for second shim | strcat('_registered_slice',num2str(slice-1))

% Acquire the data and corresponding tes
[scan,tes,GRE_info,file_name] = get_data(img_name_append, info_name_append, sequence, datafile, Pt, slice, relevant_echoes);

% Determine which dynamics should be included in the library of references
clear WFlib_partitions;
WFlib_partitions = reshape(discard_dyns+1:discard_dyns+1+num_stacks*num_avgs-1,num_stacks*num_avgs,[]);

% Compute body mask
figure;
abs_data_echo_all = squeeze(sum(abs(scan(:,:,1:end,:)),4));
abs_data_echo_all_norm = abs_data_echo_all./ max(abs_data_echo_all,[],'all');
sd_scan = std(squeeze(abs_data_echo_all_norm),0,3);

% Thresholding based on magnitude image of the first GRE
body_mask_mGRE = zeros([size(scan,1) size(scan,2)]);
body_mask_mGRE(abs(scan(:,:,1,1)) > 1.0*mean(abs(scan(:,:,1,1)),'all')) = 1;
body_mask_mGRE = imclose(body_mask_mGRE, strel('disk', 10));
subplot(1,2,1); imagesc(sd_scan, [0 0.10]); axis image off; colormap gray;
subplot(1,2,2); imagesc(body_mask_mGRE); axis image off; colormap gray;


%% Express the first dynamic into water and fat components (W, F, R2*, and dw0)

% Set water/fat recon parameters
visualize = 1;                                          % provide a visualization ('0' for no and '1' for yes)
waterfatparams.R2star_range = [0 100 21];               % [minimum | maximum | number_for_quantization]
waterfatparams.off_resonance_range = [-1500 1500 501];  % [minimum | maximum | number_for_quantization]
waterfatparams.itr = 30;                                % Number of graph cut iterations 
waterfatparams.lambda = 0.05;                          % Regularization parameter (smoothening of static off-resonance: smaller is less smoothening)
waterfatparams.threshold = 0.05;                        % 
waterfatparams.num_fat_peaks = 3;                 	    % Coose an x-peak model, options are:
                                                        % [ 1 | 3 | 4 | 5 | 6 | 7 | 9 ]

% Apply water/fat separation to the dynamics in the library of references
clear Wlib Flib dw0lib R2starlib
for REPPConditionIx = 1:size(WFlib_partitions,2)
    [Wlib(:,:,:,REPPConditionIx),Flib(:,:,:,REPPConditionIx),dw0lib(:,:,:,REPPConditionIx),R2starlib(:,:,:,REPPConditionIx),fatmodel] = perform_water_fat_separation( ...
        squeeze(scan(:,:,2:end,WFlib_partitions(:,REPPConditionIx)')),GRE_info,file_name,frequency_system,b0,prc,waterfatparams,visualize);
end


%% Masking

% Compute the fat, body and no signal mask with current parameters
[fat_mask,body_mask,brain_mask] = segmentBrainRadial(scan(:,:,1:end,WFlib_partitions(1)), body_mask_mGRE);

% Provide a visualization
figure;
subplot(1,3,1); imagesc(body_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('body mask');
subplot(1,3,2); imagesc(fat_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('fat mask');
subplot(1,3,3); imagesc(brain_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('brain mask');

    
%% Load the characteristic motion-induced susceptibility maps

% These maps are aquired using the .m file 'MotionCorrectionMapsPolarity.m', in which 
% B0 maps are calculated, and on which Singular Value Decomposition is subsequently used
dir_MotionCorrectionMaps = [datafile,'M',num2str(Pt,'%04u'),'\Unwrapped\','motion_slice',num2str(slice),'\SVDSeparate'];

coefs = 1*[2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0];
algp.B0_motion = zeros([size(scan,1),size(scan,2),8]);
algp.B0_motion(:,:,1) = coefs(1)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_sag_neg.nii']);
algp.B0_motion(:,:,2) = coefs(2)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_sag_pos.nii']);
algp.B0_motion(:,:,3) = coefs(3)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_cor_neg.nii']);
algp.B0_motion(:,:,4) = coefs(4)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_cor_pos.nii']);
algp.B0_motion(:,:,5) = coefs(6)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_FH_neg.nii']);
algp.B0_motion(:,:,6) = coefs(7)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_FH_pos.nii']);
algp.B0_motion(:,:,7) = coefs(7)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_tra_neg.nii']);
algp.B0_motion(:,:,8) = coefs(8)*pi*niftiread([dir_MotionCorrectionMaps,'\motion_map_tra_pos.nii']);

figure;
subplot(2,4,1); imagesc(algp.B0_motion(:,:,1), [-50 50]); axis off; title('LR rotation'); colorbar; axis image; axis off;
subplot(2,4,2); imagesc(algp.B0_motion(:,:,3), [-50 50]); axis off; title('AP rotation'); colorbar; axis image; axis off;
subplot(2,4,3); imagesc(algp.B0_motion(:,:,5), [-50 50]); axis off; title('FH translation'); colorbar; axis image; axis off;
subplot(2,4,4); imagesc(algp.B0_motion(:,:,7), [-50 50]); axis off; title('FH rotation'); colorbar; axis image; axis off;
subplot(2,4,5); imagesc(algp.B0_motion(:,:,2), [-50 50]); axis off; colorbar; axis image; axis off;
subplot(2,4,6); imagesc(algp.B0_motion(:,:,4), [-50 50]); axis off; colorbar; axis image; axis off;
subplot(2,4,7); imagesc(algp.B0_motion(:,:,6), [-50 50]); axis off; colorbar; axis image; axis off;
subplot(2,4,8); imagesc(algp.B0_motion(:,:,8), [-50 50]); axis off; colorbar; axis image; axis off;


%% Temporal averaging (not applied)

% Create averaging kernel for the 3rd dimension
window = 1; % odd number recommended
kernel = ones(1, 1, 1, window) / window;

% Apply convolution along 3rd dimension
scan_T_avg = convn(scan, kernel, 'same');


%% Algorithm parameters

clear scanp lib;

% Parameter structure: algorithm parameters
algp.order = 3;                   % Polynomial / spherical background phase shift order for optimization (usually 2)
algp.method = 'spher';            % 'poly' | 'spher' | 'Legendre' | 'rbf' | 'zernike'
algp.nmiter = 5;                  % # m iterations
algp.nciter = 5;                  % # c iterations
algp.stopthresh = 1e-2;           % algorithm stopping threshold
algp.max_iters = 10;              % maximum # of iterations
algp.rad2degC = rad2degC;         % Conversion between degreeC and radians
algp.suppress_info = 0;           % Whether command window information is shown per iteration or not

% Settings for fitting paramerters
algp.Acswitch = 1;                % whether drift will be fitted for or not
algp.Acmotionswitch = 1;          % whether motion will be fitted for or not
algp.positivity_constrain = 0;    % whether temperature is constrained to be positive when motion and drift are determined
initWithPreviousEstimate = 0;     % Parameter to speed up the algorithm

% Parameter structure: algorithm masking parameters
algp.body_mask = body_mask_mGRE;
algp.Nstd_ol = 0;                 % Filter out local extreme values in the fat border
algp.neighbours = 0;              % Define local neighbourhood
algp.min_nb = 0;                  % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)

% Parameter structure: scan parameters
scanp.dim = [size(scan,1) size(scan,2)];
scanp.b0 = b0;
scanp.tes = tes(2:end);
scanp.prc = prc;


%% Perform iterative optimization to separate the temperature and drift fields
clear lib c_motion cfunc mfunc Acfunc phase_motion wts phi Acinit c_motion_init minit;
disp('Optimization of multi-echo fat-suppressed MR thermometry has started..');

% Provide which dynamics are to be implemented in the optimization
show_ix = size(scan,4)-num_avgs-discard_dyns;   % number of dynamics that will be evaluated
init_ix = num_avgs+discard_dyns+1;              % first dynamic to be evaluated
step_sz = 1;                                    % step size in dynamics

% Voxel location for monitoring of temperature during optimization
px = 85;
py = 90;

tic
for ii = init_ix:step_sz:init_ix+show_ix-1
    
    % Evaluate a single dynamic
    dyn_img = scan_T_avg(:,:,2:end,ii);
    scanp.dyn_img_mask = scan_T_avg(:,:,:,ii);

    % Initialize motion parameters for the dynamic being evaluated
    for i = 1:size(WFlib_partitions,1)
        % Convert T_rigid to c_motion
        T_rigid_measured = T_rigid_temperature(:,ii)-T_rigid_temperature(:,WFlib_partitions(i,1));
        c_motion_single_refdyn = zeros([8 size(T_rigid_measured,2)]);
        if T_rigid_measured(2,:) < 0
            c_motion_single_refdyn(1,:) = -1*T_rigid_measured(2,:);
            c_motion_single_refdyn(2,:) = 0;
        else
            c_motion_single_refdyn(1,:) = 0;
            c_motion_single_refdyn(2,:) = T_rigid_measured(2,:);
        end
        if T_rigid_measured(3,:) < 0
            c_motion_single_refdyn(3,:) = -1*T_rigid_measured(3,:);
            c_motion_single_refdyn(4,:) = 0;
        else
            c_motion_single_refdyn(3,:) = 0;
            c_motion_single_refdyn(4,:) = T_rigid_measured(3,:);
        end
        if T_rigid_measured(6,:) < 0
            c_motion_single_refdyn(5,:) = -1*T_rigid_measured(6,:);
            c_motion_single_refdyn(6,:) = 0;
        else
            c_motion_single_refdyn(5,:) = 0;
            c_motion_single_refdyn(6,:) = T_rigid_measured(6,:);
        end
        if T_rigid_measured(1,:) < 0
            c_motion_single_refdyn(7,:) = -1*T_rigid_measured(1,:);
            c_motion_single_refdyn(8,:) = 0;
        else
            c_motion_single_refdyn(7,:) = 0;
            c_motion_single_refdyn(8,:) = T_rigid_measured(1,:);
        end
        scanp.T_rigid(:,i) = c_motion_single_refdyn;
    end

    % Update W/F separated reference data
    lib.Wlib = Wlib(:,:,:,1);
    lib.Flib = Flib(:,:,:,1);
    lib.dw0lib = dw0lib(:,:,:,1);
    lib.R2starlib = R2starlib(:,:,:,1);
    lib.fatmodel = fatmodel;

    disp(['Dynamic ' num2str(ii)]);

    if initWithPreviousEstimate && ii == init_ix
        % Optimization of the first dynamic
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii),phase_motion(:,:,ii),c_motion(:,ii)] = ...
            OptimizeHIMMSVD(dyn_img, algp, scanp, lib,algp.Acswitch,algp.Acmotionswitch,1,0,zeros(scanp.dim),0);
    elseif initWithPreviousEstimate && ii > init_ix
        % Optimization of remainder dynamics: set initialization of temperature
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii),phase_motion(:,:,ii),c_motion(:,ii)] = ...
            OptimizeHIMMSVD(dyn_img, algp, scanp, lib,algp.Acswitch,algp.Acmotionswitch,1,0,mfunc(:,:,end),0);
    else
        % Phi switch  off, m switch on (solve for temperature without Tx/Rx gain)
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii),phase_motion(:,:,ii),c_motion(:,ii)] = ...
            OptimizeHIMMSVD(dyn_img, algp, scanp, lib,algp.Acswitch,algp.Acmotionswitch,1,zeros(scanp.dim),0);
    end

    % Give visual update during optimization
    close all; 
    figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4]);
    subplot(1,3,1); 
    imagesc(mfunc(:,:,ii)*algp.rad2degC.*brain_mask, [-1 1]); hold on;
    plot(px,py,'marker','x','MarkerSize',20,'LineWidth',5,'Color',[0, 0.4470, 0.7410]);
    axis image; axis off; colorbar; title(ii);
    subplot(1,3,2); 
    current_dynamic = Acfunc(:,:,end); drift_voxels = current_dynamic(body_mask == 1);
    vmax_current = max(drift_voxels)+1; vmin_current = min(drift_voxels)-1;
    imagesc(Acfunc(:,:,ii).*body_mask, [vmin_current vmax_current]);
    axis image; axis off; colorbar; title(ii);
    subplot(1,3,3); 
    plot(1:ii,squeeze(mfunc(py,px,1:ii)*algp.rad2degC),'Color',[0, 0.4470, 0.7410]); 
    axis square; ylim([-0.5 0.5]); xlim([1 size(scan,4)]);
    drawnow;
    
end
toc


%% Quick visualization

figure;
p1y = 96; p1x = 75;
p2y = 137; p2x = 109;
subplot(1,2,1); imagesc(mean(mfunc(:,:,end-10:end),3)*rad2degC.*brain_mask, [-1 1]); axis image; axis off; colorbar; hold on
plot(p1x,p1y,'marker','x','MarkerSize',20,'LineWidth',5,'Color',[0, 0.4470, 0.7410]);
plot(p2x,p2y,'marker','x','MarkerSize',20,'LineWidth',5,'Color',[0.8500, 0.3250, 0.0980]);
subplot(1,2,2); plot(squeeze(mfunc(p1y,p1x,:))*rad2degC, LineWidth=2); hold on;
plot(squeeze(mfunc(p2y,p2x,:))*rad2degC, LineWidth=2); hold on;
ylim([-1 1]); xlim([1 size(mfunc,3)]); axis square; 