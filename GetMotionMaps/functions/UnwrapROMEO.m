function [unwrapped_phase_scan,B0_scan] = UnwrapROMEO(scan,body_mask,tes,GRE_info,template)
% Phase unwrapping using ROMEO to prevent having to deal with wrap-around
% artifacts. Gitub to ROMEO: https://github.com/korbinian90/ROMEO

% Set parameter settings and define the output path
parameters.output_dir = fullfile(['C:\Users\mkikken\OneDrive - UMC Utrecht\Documenten\PhD\ROMEO_tmp']);
parameters.mask = 'nomask';
if ~isfolder(parameters.output_dir)
    mkdir(parameters.output_dir);
end

% Perform the unwrapping procedure as a temporal series
% Load the magnitude and phase images with the echo in the fourth dimension: [nz nx ny nechoes]
parameters.mag = permute(abs(scan(:,:,:,:)).*repmat(body_mask,[1,1,size(scan,3),size(scan,4)]), [4 1 2 3]);
phase = permute(angle(scan(:,:,:,:)).*repmat(body_mask,[1,1,size(scan,3),size(scan,4)]), [4 1 2 3]);

% Load a body mask [nz nx ny]
parameters.mask = permute(repmat(int8(body_mask),[1,1,size(scan,4)]),[3 1 2]);

% Define parameters
parameters.calculate_B0 = true; 
parameters.phase_offset_correction = 'off';
%parameters.phase_offset_smoothing_sigma_mm = [50,50,50];
parameters.additional_flags = '-i';
parameters.TE = tes'*1e3; % echo times in [ms]
parameters.voxel_size = [GRE_info.imgdef.slice_thickness.uniq GRE_info.imgdef.pixel_spacing_x_y.uniq];

% Unwrap the image
[unwrapped, B0] = ROMEO(phase, parameters);
unwrapped(isnan(unwrapped)) = 0;

% Convert back to [nx ny nechoes] and put into array containing all dynamics
unwrapped_phase_scan = squeeze(permute(unwrapped, [2 3 4 1]));
B0_scan = squeeze(permute(B0, [2 3 1]));


% % Perform the unwrapping procedure dynamic-by-dynamic
% unwrapped_phase_scan = zeros(size(scan)); B0_scan = zeros(size(squeeze(scan(:,:,1,:))));
% num_dynamics = size(scan,4);
% for dyn = 1:num_dynamics
% 
%     % Provide a number in the command window such that the user knows the progress
%     disp(['   Dynamic ' num2str(dyn) ' of ' num2str(num_dynamics)]);
% 
%     % Load the magnitude and phase images with the echo in the fourth dimension: [nz nx ny nechoes]
%     parameters.mag = permute(abs(scan(:,:,:,dyn)).*repmat(body_mask,[1,1,size(scan,3)]), [4 1 2 3]);
%     phase = permute(angle(scan(:,:,:,dyn)).*repmat(body_mask,[1,1,size(scan,3)]), [4 1 2 3]);
% 
%     % Provide template (temporal average)
%     parameters.template = template;
% 
%     % Load a body mask [nz nx ny]: based on the water/fat separated image
%     parameters.mask = permute(int8(body_mask), [3 1 2]);
% 
%     % Define parameters
%     parameters.calculate_B0 = true; 
%     parameters.phase_offset_correction = 'off'; % 'off' | 'bipolar'
%     parameters.additional_flags = '-i'; %  -g -i
%     %parameters.weights = 'romeo';
%     %parameters.B0_phase_weighting = 'TEs';
%     parameters.TE = tes'*1e3; % echo times in [ms]
%     parameters.voxel_size = [GRE_info.imgdef.slice_thickness.uniq GRE_info.imgdef.pixel_spacing_x_y.uniq];
% 
%     % Unwrap the image
%     [unwrapped, B0] = ROMEO(phase, parameters);
%     unwrapped(isnan(unwrapped)) = 0;
% 
%     % Convert back to [nx ny nechoes] and put into array containing all dynamics
%     unwrapped_phase_scan(:,:,:,dyn) = squeeze(permute(unwrapped, [2 3 4 1]));
%     B0_scan(:,:,dyn) = squeeze(permute(B0, [2 3 1]));
% end
% 
% rmdir(parameters.output_dir, 's') % remove the temporary ROMEO output folder
% disp(['Done!']);

end