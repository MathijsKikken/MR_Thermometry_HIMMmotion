function [Ac,mask_outliers] = NH2DReconPreparation(wts,c_motion,lib,scanp,algp,fat_idx,Acmotionswitch)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% NH2DReconPreparation: perpare the drift field in the fat layer to make it
% appropriate for near harmonic 2D reconstruction
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   wts:                baseline weights
%   c_motion:           motion fields coefficients fit
%   lib:                baseline library structure
%   scanp:              scan params structure
%   algp:               algorithm params structure
%   fat_idx:            fat mask
%   Acmotionswitch      integer specifying whether we should include a fit
%                       for motion in the algorithm
%
% Output: 
%   Ac:                 frequency shift due to drift, estimated using
%                       phase-subtraction and corrected for motion
%   mask_outliers:      mask out regions that should not be taken into
%                       account for near-harmonic 2D reconstruction

    % Make sure the reference is based on the weighted combination of library baseline images
    LibReference = zeros(size(lib.imgs_ref(:,:,:,1)));
    for lib_ix = 1:size(lib.imgs_ref,4)
        LibReference = LibReference + wts(lib_ix)*lib.imgs_ref(:,:,:,lib_ix);
    end             
    
    % Add the motion contributions onto the phase accumulation map
    scan_motion_corr = zeros(scanp.dim);
    motion_contribution = zeros(scanp.dim);
    if Acmotionswitch
        for i = 1:size(c_motion,1)
            motion_contribution = motion_contribution + algp.B0_motion(:,:,i)*c_motion(i);
        end
    end
    for e = 1:size(scanp.tes,1)
        corrected_phase = angle(lib.imgs(:,:,e)) - motion_contribution*scanp.tes(e);
        scan_motion_corr(:,:,e) = abs(lib.imgs(:,:,e)).*exp(1j*corrected_phase);
    end
    phase_acc_motion_corrected = angle(scan_motion_corr.*conj(LibReference));
    
    % Average over echoes
    avg_map_TEw = zeros(size(phase_acc_motion_corrected(:,:,1)));
    for echo = 1:size(phase_acc_motion_corrected,3)
       temp = squeeze(phase_acc_motion_corrected(:,:,echo))/scanp.tes(echo);
       avg_map_TEw = avg_map_TEw + temp;
    end
    avg_map_TEw = avg_map_TEw / (size(phase_acc_motion_corrected,3)-1+1);
    Ac = avg_map_TEw.*fat_idx;
    
    % Do not take regions into account where B0 was not determined for motion correction
    mask_outliers = fat_idx;
    mask_B0 = ones(size(fat_idx));
    mask_B0(sum(algp.B0_motion,3) == 0) = 0;
    mask_outliers = mask_outliers.*mask_B0;
    fat_idx = fat_idx.*mask_outliers;
    
    % Perform median filter, where values outside of the fat/body mask are excluded
    size_kernel = algp.median_filter_size_NHrecon;
    Ac_premedfilter = Ac;
    Ac_premedfilter(fat_idx == 0) = nan;
    Ac = nanmedfilt2(Ac_premedfilter,size_kernel);
    Ac(isnan(Ac)) = 0;

end