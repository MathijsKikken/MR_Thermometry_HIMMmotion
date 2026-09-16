%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% 
%  HIMM: Harmonic Initialized Model-based Multi-echo  
%  SVD: Singular Value Decomposition (used as a motion correction scheme)
%  Sub script to extract the characteristic motion fields
%
%  This algorithm code determines the B0 field of all individual dynamics.
%  Individual B0 fields are extracted from multi-echo data through the
%  ROMEO unwrapping algorithm (https://github.com/korbinian90/ROMEO)
%
%  After calculation of the individual B0 fields, every dynamic is compared
%  with a near-neutral reference dynamic from the same scan. Positive and
%  negative motion states are separated. A separate singular value
%  decomposition is then performed for each motion type and direction to
%  obtain direction-specific characteristic motion-induced B0 field maps.
%
%  Complete main script, including select_signed_motion and esa below.
%  Requires your existing get_data, UnwrapROMEO and SVD_motion functions,
%  ROMEO, Image Processing Toolbox and Statistics and Machine Learning Toolbox.
%  Motion parameters are assumed to be in degrees (rotations) and mm
%  (translations). The reference is closest to zero in the selected parameter;
%  verify that this corresponds to the intended neutral head position.
%  Original R2 weighting of reconstructed components is retained.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
%
% Do not reproaduce, distribute, or modify without proper citation according
% to license file

clear all; close all; clc;

% Add the software (ROMEO phase unwrapping and B0 calculation algorithm) to the directory
addpath(genpath('.\Code\GetMotionMaps\functions'));
addpath(genpath('.\Code\MRITools')); % ROMEO


%% Provide subject data

% Patient definition
datafile  = '.\Data\';
Pt        = 1;                % subject number
slice     = 1;                % usually set to 1, except for M2D acquisitions

% Specify which echoes are used for the reconstruction, and the appendix of
% the image and info file (after '_img' or _info', e.g. '_img_shim2')
relevant_echoes = [1,2:5];
global_append = '_registered_slice0'; % Prior to estimation of motion maps, all maps are registered to the first dynamics in the transverse plane
                                      % '' for raw data
                                      % '_registered_slice0' for data registed to the first dynamic


%% Load the rigid transformation parameters of the dynamic and motion scan(s) 
% Transformation parameters of the motion maps are used to define a library
% of common motion states. The transformation parameters of the dynamic scan
% is subsequently used to initialize the spatial motion map during optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Rotation in the sagittal plane (nodding)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the data (incl. T_rigid) in which motion was deliberately induced
img_name_append = strcat('_motion_sag',global_append); info_name_append = '_motion_sag';
[scan_sag,tes_sag,GRE_info_sag,~] = get_data(img_name_append, info_name_append, 'baseline', datafile, Pt, slice, relevant_echoes);
T_rigid_loaded_sag = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_motion_sag','.nii.gz']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Rotation in the coronal plane (left/right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the data (incl. T_rigid) in which motion was deliberately induced
img_name_append = strcat('_motion_cor',global_append); info_name_append = '_motion_cor';
[scan_cor,tes_cor,GRE_info_cor,~] = get_data(img_name_append, info_name_append, 'baseline', datafile, Pt, slice, relevant_echoes);
T_rigid_loaded_cor = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_motion_cor','.nii.gz']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Translation in the feet/head (FH) direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the data (incl. T_rigid) in which motion was deliberately induced
img_name_append = strcat('_motion_FH',global_append); info_name_append = '_motion_FH';
[scan_FH,tes_FH,GRE_info_FH,~] = get_data(img_name_append, info_name_append, 'baseline', datafile, Pt, slice, relevant_echoes);
T_rigid_loaded_FH = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_motion_FH','.nii.gz']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Rotation in the transverse plane (left/right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the data (incl. T_rigid) in which motion was deliberately induced
img_name_append = strcat('_motion_tra',global_append); info_name_append = '_motion_tra';
[scan_tra,tes_tra,GRE_info_tra,~] = get_data(img_name_append, info_name_append, 'baseline', datafile, Pt, slice, relevant_echoes);
T_rigid_loaded_tra = niftiread([datafile 'M',num2str(Pt,'%04u'),'\T_rigid_motion_tra','.nii.gz']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Get a mask of the brain using thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thresholding based on magnitude image of the first GRE
mask = zeros([size(scan_sag,1) size(scan_sag,2)]);
mask(abs(scan_sag(:,:,1,1)) > 1.0*mean(abs(scan_sag(:,:,1,1)),'all')) = 1;
mask = imclose(mask, strel('disk', 10));
figure; imagesc(mask); axis image; axis off; title('mask');


%% Calculate the B0 map for every motion state individually (using ROMEO)

% Perform phase unwrapping using ROMEO - check if phase unwrapping was already performed
output_dir = [datafile,'\M',num2str(Pt,'%04u'),'\Unwrapped\','motion_slice',num2str(slice)];
if ~exist(output_dir, 'dir')

    % Phase unwrapping using ROMEO
    disp('Unwrapping data with Rot. around LR axis');
    [unwrapped_phase_scan_sag,B0_sag] = UnwrapROMEO(scan_sag(:,:,2:end,:),mask,tes_sag(2:end,:),GRE_info_sag,1);
    disp('Unwrapping data with Rot. around AP axis');
    [unwrapped_phase_scan_cor,B0_cor] = UnwrapROMEO(scan_cor(:,:,2:end,:),mask,tes_cor(2:end,:),GRE_info_cor,1);
    disp('Unwrapping data with FH translations');
    [unwrapped_phase_scan_FH,B0_FH] = UnwrapROMEO(scan_FH(:,:,2:end,:),mask,tes_FH(2:end,:),GRE_info_FH,1);
    disp('Unwrapping data with Rot. around FH axis');
    [unwrapped_phase_scan_tra,B0_tra] = UnwrapROMEO(scan_tra(:,:,2:end,:),mask,tes_tra(2:end,:),GRE_info_tra,1);
    
    % Save for future use
    mkdir(output_dir);
    %niftiwrite(unwrapped_phase_scan,[output_dir,'\unwrapped_phase_scan.nii']);
    niftiwrite(B0_sag,[output_dir,'\B0_sag.nii']);
    niftiwrite(B0_cor,[output_dir,'\B0_cor.nii']);
    niftiwrite(B0_FH,[output_dir,'\B0_FH.nii']);
    niftiwrite(B0_tra,[output_dir,'\B0_tra.nii']);
else
    % Phase unwrapping was already performed, so we load from the NiFTI files
    %niftiread([output_dir,'\unwrapped_phase_scan.nii']);
    B0_sag = niftiread([output_dir,'\B0_sag.nii']);
    B0_cor = niftiread([output_dir,'\B0_cor.nii']);
    B0_FH = niftiread([output_dir,'\B0_FH.nii']);
    B0_tra = niftiread([output_dir,'\B0_tra.nii']);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate B0 changes relative to a neutral state and keep positive and
% negative motion separate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rows of T_rigid containing the deliberately induced motion:
%   sagittal scan: rotation around LR axis
%   coronal scan:  rotation around AP axis
%   transverse scan: rotation around FH axis
%   FH scan: translation in FH direction
motion_rows = [2, 3, 1, 6];

% For each scan, select the dynamic closest to zero in the relevant motion
% parameter as its reference. Each other dynamic is compared only with this
% reference. The sign therefore describes motion away from the neutral state,
% rather than the arbitrary ordering of a pair returned by nchoosek.
[B0_sag_rel,T_rigid_sag_rel,motion_sign_sag,ref_sag] = ...
    select_signed_motion(B0_sag,T_rigid_loaded_sag(1:6,:),motion_rows(1));
[B0_cor_rel,T_rigid_cor_rel,motion_sign_cor,ref_cor] = ...
    select_signed_motion(B0_cor,T_rigid_loaded_cor(1:6,:),motion_rows(2));
[B0_tra_rel,T_rigid_tra_rel,motion_sign_tra,ref_tra] = ...
    select_signed_motion(B0_tra,T_rigid_loaded_tra(1:6,:),motion_rows(3));
[B0_FH_rel,T_rigid_FH_rel,motion_sign_FH,ref_FH] = ...
    select_signed_motion(B0_FH,T_rigid_loaded_FH(1:6,:),motion_rows(4));

fprintf('Reference dynamics: sagittal %d, coronal %d, transverse %d, FH %d\n', ...
    ref_sag,ref_cor,ref_tra,ref_FH);

% Combine only the motion parameters for the overview plot below.
% The B0 maps themselves remain separated.
T_rigid = cat(2,T_rigid_sag_rel,T_rigid_cor_rel,T_rigid_tra_rel,T_rigid_FH_rel);

% Create the eight SVD inputs directly from the four relative B0 series.
% No combined B0 array is used as input to SVD_motion.
B0_groups = { ...
    B0_sag_rel(:,:,motion_sign_sag > 0), ...
    B0_sag_rel(:,:,motion_sign_sag < 0), ...
    B0_cor_rel(:,:,motion_sign_cor > 0), ...
    B0_cor_rel(:,:,motion_sign_cor < 0), ...
    B0_tra_rel(:,:,motion_sign_tra > 0), ...
    B0_tra_rel(:,:,motion_sign_tra < 0), ...
    B0_FH_rel(:,:,motion_sign_FH > 0), ...
    B0_FH_rel(:,:,motion_sign_FH < 0)};

T_groups = { ...
    T_rigid_sag_rel(:,motion_sign_sag > 0), ...
    T_rigid_sag_rel(:,motion_sign_sag < 0), ...
    T_rigid_cor_rel(:,motion_sign_cor > 0), ...
    T_rigid_cor_rel(:,motion_sign_cor < 0), ...
    T_rigid_tra_rel(:,motion_sign_tra > 0), ...
    T_rigid_tra_rel(:,motion_sign_tra < 0), ...
    T_rigid_FH_rel(:,motion_sign_FH > 0), ...
    T_rigid_FH_rel(:,motion_sign_FH < 0)};


%% Plot rotation/translation parameters 
figure;
subplot(4,1,1);
plot(T_rigid(2,:),linewidth=2); title('LR axis rotation')
subplot(4,1,2);
plot(T_rigid(3,:),linewidth=2,color='r'); title('AP axis rotation')
subplot(4,1,3);
plot(T_rigid(6,:),linewidth=2,color='black'); title('FH translation')
subplot(4,1,4);
plot(T_rigid(1,:),linewidth=2,color='g'); title('FH axis rotation')


%% Create output folders for the separate SVD analyses

% Output folder: where figures are stored
output_dir_separate = [output_dir,'\SVDSeparate'];
if ~exist(output_dir_separate, 'dir')
    mkdir(output_dir_separate);
end
output_dir_figs = [output_dir_separate,'\Figures'];
if ~exist(output_dir_figs, 'dir')
    mkdir(output_dir_figs);
end

% Plot motion parameters
TrigidPlot = {T_rigid_loaded_sag,T_rigid_loaded_cor,T_rigid_loaded_tra,T_rigid_loaded_FH};
for i = 1:4
    T_plot = TrigidPlot{i};
    figure; subplot(1,1,1);
    plot(T_plot(2,:),linewidth=4); hold on;
    plot(T_plot(3,:),linewidth=4,color='r'); hold on;
    plot(T_plot(6,:),linewidth=4,color='black'); hold on;
    Trigid_plot = plot(T_plot(1,:),linewidth=4,color='g');
    max_yval = max(abs(T_plot([1,2,3,6],:)),[],'all')+0.1;
    xlim([1 max(2,size(T_plot,2))]); ylim([-max_yval max_yval]);
    x0=10; y0=10; width=1000;height=400;
    set(gcf,'position',[x0,y0,width,height]); hold off;
    %saveas(Trigid_plot,[output_dir_figs,'\Trigid_plot_motion',num2str(i),'.svg']);
end
close all;


%% Perform a separate SVD for every motion direction and sign

% Define the eight sub-scans: one motion type and one sign per analysis.
motion_types = [2,-2,3,-3,1,-1,6,-6];
motion_rows_fit = abs(motion_types);
motion_labels = {'sag_pos','sag_neg','cor_pos','cor_neg', ...
                 'tra_pos','tra_neg','FH_pos','FH_neg'};
plot_titles = {'Rot. LR +','Rot. LR -','Rot. AP +','Rot. AP -', ...
               'Rot. FH +','Rot. FH -','Trans. FH +','Trans. FH -'};
colors = ["#0072BD";"#56B4E9";"#D95319";"#EDB120"; ...
          "#008000";"#77AC30";"#000000";"#7E2F8E"];

num_vectors_max = 5;  % Maximum number of LSVs evaluated per sub-scan
min_R2 = 0.01;        % Minimum R2 for including an LSV
fit_bias = 1;         % 0: fit through origin; 1: include an intercept

% Unavailable fits remain NaN rather than appearing as zero field changes.
B0_motion = nan(size(B0_sag_rel,1),size(B0_sag_rel,2),numel(motion_types));

for motion_ix = 1:numel(motion_types)
    % Select the explicitly separated array for this SVD.
    B0_sub = B0_groups{motion_ix};
    T_sub = T_groups{motion_ix};
    motion_row = motion_rows_fit(motion_ix);
    x_data = T_sub(motion_row,:).';

    if numel(x_data) < 2 + fit_bias || numel(unique(x_data)) < 2
        warning('Skipping %s: insufficient motion states or motion variation.',motion_labels{motion_ix});
        continue
    end
    assert(all(isfinite(B0_sub(:))) && all(isfinite(x_data)), ...
        'Nonfinite B0 or motion values in %s.',motion_labels{motion_ix});
    B0_motion(:,:,motion_ix) = 0;

    % The SVD basis is now specific to this motion type and direction.
    [U,S,V] = SVD_motion(B0_sub);
    num_vectors = min([num_vectors_max,size(U,2),size(S,1),size(S,2),size(V,2)]);
    B0_size_sub = size(B0_sub);

    % S*V' contains the coefficient of each LSV for every B0 map.
    weights = S(1:num_vectors,1:num_vectors) * ...
              V(:,1:num_vectors).';

    % Convert the retained left singular vectors back into spatial maps.
    LSVs = zeros(B0_size_sub(1),B0_size_sub(2),num_vectors);
    for vector_ix = 1:num_vectors
        LSVs(:,:,vector_ix) = reshape(U(:,vector_ix), ...
            [B0_size_sub(1),B0_size_sub(2)]);

        figure;
        fig_LSV = imagesc(LSVs(:,:,vector_ix),[-0.015 0.015]);
        axis image; axis off; colormap(esa);
        saveas(fig_LSV,[output_dir_figs,'\LSV_',motion_labels{motion_ix}, ...
            '_',num2str(vector_ix),'.svg']);
        close(gcf);
    end

    % Correlate the sub-scan-specific LSV weights with its motion parameter.
    biass = zeros(1,num_vectors);
    dcs = zeros(1,num_vectors);
    R2s = zeros(1,num_vectors);
    for vector_ix = 1:num_vectors
        if fit_bias == 1
            mdl = fitlm(x_data,weights(vector_ix,:).');
            biass(vector_ix) = mdl.Coefficients.Estimate(1);
            dcs(vector_ix) = mdl.Coefficients.Estimate(2);
        else
            mdl = fitlm(x_data,weights(vector_ix,:).','Intercept',false);
            dcs(vector_ix) = mdl.Coefficients.Estimate(1);
        end
        R2s(vector_ix) = mdl.Rsquared.Ordinary;
    end

    % Visualize up to the first five fits and the last evaluated vector.
    vectors_to_plot = unique([1:min(5,num_vectors),num_vectors],'stable');
    figure('units','normalized','outerposition',[0 0 1 1]);
    ylim_vals = max(abs(weights),[],'all');
    if ylim_vals == 0
        ylim_vals = 1;
    end
    x_axis_motion = linspace(min(x_data),max(x_data),100);
    for plot_ix = 1:numel(vectors_to_plot)
        vector_plot = vectors_to_plot(plot_ix);
        subplot(1,numel(vectors_to_plot),plot_ix);
        p = plot(x_data,weights(vector_plot,:).','o','MarkerSize',3, ...
            'Color',colors(motion_ix,:)); hold on;
        p.MarkerFaceColor = colors(motion_ix,:);
        fig_fit = plot(x_axis_motion,biass(vector_plot) + ...
            x_axis_motion*dcs(vector_plot),'Color',"#EDB120",'linewidth',1);
        title(sprintf('LSV %d, R^2 = %.2f',vector_plot,R2s(vector_plot)));
        ylim([-ylim_vals ylim_vals]); axis square;
    end
    saveas(fig_fit,[output_dir_figs,'\LinearFit_',motion_labels{motion_ix},'.svg']);
    close(gcf);

    % Reconstruct the characteristic field at +1 or -1 degree/mm.
    target_motion = sign(motion_types(motion_ix));
    for vector_ix = 1:num_vectors
        target_weight = target_motion*dcs(vector_ix);
        B0_motion(:,:,motion_ix) = B0_motion(:,:,motion_ix) + ...
            LSVs(:,:,vector_ix)*R2s(vector_ix)*target_weight;
    end

    niftiwrite(B0_motion(:,:,motion_ix), ...
        [output_dir_separate,'\motion_map_',motion_labels{motion_ix},'.nii']);
end

% Show the eight independently estimated characteristic motion maps.
figure('units','normalized','outerposition',[0 0 1 1]);
for motion_ix = 1:numel(motion_types)
    subplot(2,4,motion_ix);
    fig_char_motion = imagesc(B0_motion(:,:,motion_ix),[-10 10]);
    axis image; axis off; title(plot_titles{motion_ix}); colormap(esa);
end
saveas(fig_char_motion,[output_dir_figs,'\CharacteristicMotionMaps.svg']);


%% Select signed motion relative to the neutral dynamic
function [B0_rel,T_rel,motion_sign,ref_idx] = select_signed_motion(B0,T_rigid,motion_row)

assert(size(B0,3) == size(T_rigid,2), ...
    'B0 dynamics and rigid-motion parameter counts must match.');
assert(all(isfinite(T_rigid(:))), ...
    'Motion parameters must be finite.');

% Select the acquired dynamic closest to the median motion parameter.
motion_values = T_rigid(motion_row,:);
[~,ref_idx] = min(abs(motion_values - median(motion_values)));

% Calculate differences relative to this reference.
use_idx = setdiff(1:size(B0,3),ref_idx,'stable');

B0_rel = B0(:,:,use_idx) - B0(:,:,ref_idx);
T_rel = T_rigid(:,use_idx) - T_rigid(:,ref_idx);

% Classify motion relative to the selected reference.
motion_sign = sign(T_rel(motion_row,:)).';

% Remove states with zero displacement in the relevant parameter.
keep = motion_sign ~= 0;

B0_rel = B0_rel(:,:,keep);
T_rel = T_rel(:,keep);
motion_sign = motion_sign(keep);

end


%% colormap function
function h = esa(m)
%ESA    Color map providing ESA corporate design colours.
%       ESA(M) returns an M-by-3 matrix of the colormap.
if nargin < 1
    m = size(get(gcf,'colormap'),1); 
else
    if isempty(m), m = size(get(gcf,'colormap'),1); end
    if ischar(m), m = str2double(m); end
end
h = [     
     0.356860000000000011     0.466669999999999974     0.800000000000000044;
     0.387490000000000001     0.492070000000000007     0.809520000000000017;
     0.418109999999999982     0.517460000000000031     0.819049999999999945;
     0.448740000000000028     0.542860000000000009     0.828570000000000029;
     0.479360000000000008     0.568259999999999987     0.838099999999999956;
     0.509990000000000054     0.593650000000000011      0.84762000000000004;
     0.540610000000000035     0.619049999999999989     0.857140000000000013;
      0.57123999999999997     0.644449999999999967     0.866670000000000051;
     0.601870000000000016     0.669839999999999991     0.876190000000000024;
     0.632489999999999997     0.695239999999999969     0.885709999999999997;
     0.663120000000000043     0.720639999999999947     0.895240000000000036;
     0.693740000000000023     0.746029999999999971     0.904760000000000009;
     0.724369999999999958     0.771429999999999949     0.914290000000000047;
      0.75499000000000005     0.796830000000000038      0.92381000000000002;
     0.785619999999999985     0.822219999999999951     0.933329999999999993;
     0.816250000000000031      0.84762000000000004     0.942860000000000031;
     0.846870000000000012     0.873020000000000018     0.952380000000000004;
     0.877499999999999947     0.898410000000000042     0.961899999999999977;
     0.908120000000000038      0.92381000000000002     0.971430000000000016;
     0.938749999999999973     0.949209999999999998     0.980949999999999989;
     0.969369999999999954     0.974600000000000022     0.990480000000000027;
                        1                        1                        1;
                        1     0.988589999999999969     0.961679999999999979;
                        1     0.977180000000000049     0.923359999999999959;
                        1     0.965770000000000017     0.885050000000000003;
                        1     0.954359999999999986     0.846729999999999983;
                        1     0.942949999999999955     0.808409999999999962;
                        1     0.931549999999999989     0.770090000000000052;
                        1     0.920139999999999958     0.731770000000000032;
                        1     0.908730000000000038     0.693450000000000011;
                        1     0.897320000000000007     0.655139999999999945;
                        1     0.885909999999999975     0.616820000000000035;
                        1     0.874500000000000055     0.578500000000000014;
                        1     0.863090000000000024     0.540179999999999993;
                        1     0.851679999999999993     0.501859999999999973;
                        1     0.840269999999999961     0.463550000000000018;
                        1     0.828860000000000041     0.425229999999999997;
                        1      0.81745000000000001     0.386909999999999976;
                        1     0.806050000000000044     0.348590000000000011;
                        1     0.794640000000000013      0.31026999999999999;
                        1     0.783229999999999982     0.271950000000000025;
                        1     0.771819999999999951     0.233639999999999987;
                        1      0.76041000000000003     0.195319999999999994;
];
end
