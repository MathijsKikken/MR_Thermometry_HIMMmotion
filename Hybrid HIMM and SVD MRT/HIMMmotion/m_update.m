function m = m_update(imgsw,Zw,m,tes,prc,mask)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% m_update: function to fit PRF shift due to heating
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgsw:              Dynamic multi-echo images water component (Nx x Ny x Necho)
%   Zw:                 Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   m:                  initial guess
%   tes:                echo times
%   prc:                direction of precession
%   body_mask:          region consisting of brain and skull
%   signal_mask:        region without signal
%   mask:               masked heat
%
% Output: 
%   m:                  dim x dim fit for temperature

% Mask measured and modeled data
brain_mask_TEs = repmat(mask,1,1,length(tes));
mask_list = brain_mask_TEs == 1;
imgsw_masked = imgsw;
imgsw_masked = imgsw_masked.*mask_list;
Zw_masked = Zw;
Zw_masked = Zw_masked.*mask_list;

% reshape images and model into vectors
innprod = Zw_masked.*conj(imgsw_masked);

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)
    hesssum = hesssum + t4(:,:,ii)./t2(:,:,ii)*tes(ii)^2;
    gradsum = gradsum - (1-2*prc)*t4(:,:,ii)*tes(ii);
end

% update m
dm = -gradsum./hesssum;
dm(isnan(dm)) = 0;
m = m + dm;

% Enforce limits
m(mask == 0) = 0;