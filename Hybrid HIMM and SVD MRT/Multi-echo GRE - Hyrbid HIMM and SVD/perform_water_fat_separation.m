function [Wlib,Flib,dw0lib,R2starlib,fatmodel] = perform_water_fat_separation(imgs,GRE_info,file_name,frequency_system,b0,prc,waterfatparams,visualize)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% perform_water_fat_separation: apply water/fat separation to baseline data
% using the robust water/fat separation algorithm using graph cuts
% (Hernando et al, MRM 2009). This algorithm can be accessed from the 
% ISMRM fat/water toolbox (https://www.ismrm.org/workshops/FatWater12/data.htm)
%
% Creators: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs: 
%   imgs:               Dynamic multi-echo images (Nx x Ny x Necho)
%   GRE_info:           info file
%   file_name:          name of the dynamic data being processed
%   frequency_system:   rate of precession
%   b0:                 b0 field strength
%   prc:                direction of precession
%   waterfatparams:     fat/water reconstruction parameters
%   visualize:          boolean whether a visualization is provided
%
% Output: 
%   Wlib:               water component of the Dixon reconstruction
%   Flib:               fat component of the Dixon reconstruction
%   dw0lib:             static off-resonance component of the Dixon reconstruction
%   R2starlib:          transverse relaxation component of the Dixon reconstruction


% Selection of fat model and construction of algorithm parameters
fatmodel = get_fatmodel(waterfatparams.num_fat_peaks);
algoParams = getReconParams(waterfatparams.R2star_range,waterfatparams.off_resonance_range,waterfatparams.itr,waterfatparams.lambda,waterfatparams.threshold,fatmodel);

disp('Separating baseline WF - using ISMRM toolbox');

% WF separation is done for each dynamic individually with image size of [ nx | ny | 1 | ncoils | nTE ]
for ii = 1:size(imgs,4)
    % Image data parameters  
    imDataParams = getImDataParams(GRE_info,imgs(:,:,:,ii),file_name,frequency_system,b0,prc);
    %algoParams.fieldmap = dw0lib(:,:,ii)/(2*pi);
    %algoParams.r2starmap = R2starlib(:,:,ii);

    % Perform Mixed-Magnitude method initialized with the Graph Cut method
    disp(['.. Executing baseline library dynamic ' num2str(ii) ' out of ' num2str(size(imgs,4))]);
    outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams,algoParams);
    %outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit(imDataParams,algoParams);
    
    % Combine data (water, fat, off-resonance and R2*) of individual dynamics
    Wlib(:,:,ii) = outParams.species(1).amps;
    Flib(:,:,ii) = outParams.species(2).amps;
    dw0lib(:,:,ii) = outParams.fieldmap*2*pi;
    R2starlib(:,:,ii) = outParams.r2starmap;

    % Smooth static off resonance map
    %dw0_presmoothing = dw0lib(:,:,ii);
    %dw0_presmoothing(body_mask_mGRE == 0) = nan;
    %dw0_postsmoothing = smoothdata2(dw0_presmoothing,"movmedian",median_filter_size,"omitnan");
    %dw0_postsmoothing(isnan(dw0_presmoothing)) = 0;
    %dw0lib(:,:,ii) = dw0_postsmoothing;
end

% Provide a visualization (magnitude and phase) of acquired water/fat separation
if visualize
    for i = 1:size(imgs,4)
        subplot(6,size(imgs,4),i);
        imagesc(abs(Wlib(:,:,i))); axis off; title(['Water | dyn ' num2str(i)]); axis image;
        subplot(6,size(imgs,4),size(imgs,4)+i); 
        imagesc(abs(Flib(:,:,i))); axis off; title(['Fat | dyn ' num2str(i)]); axis image;
        subplot(6,size(imgs,4),2*size(imgs,4)+i);
        imagesc(angle(Wlib(:,:,i))); axis off; title(['Water | dyn ' num2str(i)]); axis image;
        subplot(6,size(imgs,4),3*size(imgs,4)+i); 
        imagesc(angle(Flib(:,:,i))); axis off; title(['Fat | dyn ' num2str(i)]); axis image;
        subplot(6,size(imgs,4),4*size(imgs,4)+i); 
        imagesc(dw0lib(:,:,i)); axis off; title(['\Delta B_0 | dyn ' num2str(i)]); axis image;
        subplot(6,size(imgs,4),5*size(imgs,4)+i); 
        imagesc(R2starlib(:,:,i)); axis off; title(['R_2^* | dyn ' num2str(i)]); axis image;
    end
end

end

