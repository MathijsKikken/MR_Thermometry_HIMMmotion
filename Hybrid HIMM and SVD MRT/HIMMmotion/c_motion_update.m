function c = c_motion_update(imgs,Z,A,c,tes,prc,mask,lambda)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% c_update: function to find the coefficients on the motion fiends
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% NEWTON's METHOD GRADIENT DESCENT
%
% Inputs:
%   imgs:       Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:          Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   A:          basis function of motion fields
%   c:          previous coefficient estimate
%   tes:        echo times
%   prc:        direction of precession
%   lambda:     step size coefficient (usually 1) 
%
% Output: 
%   c:          polynomial/spherical harmonics coefficents

% Mask measured and modeled data
brain_mask_TEs = repmat(mask,1,1,length(tes));
mask_list = brain_mask_TEs == 1;
imgs_masked = imgs;
imgs_masked = imgs_masked.*mask_list;
Z_masked = Z;
Z_masked = Z_masked.*mask_list;

% reshape images and model into vectors
innprod = Z_masked.*conj(imgs_masked);
innprod = permute(innprod,[3 1 2]); 
innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
t1 = (abs(innprod));
t2 = (angle(innprod)); %wrapped w/ heat

t3 = sin(t2);
t4 = t1.*t3;

hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
    hesssum = hesssum + t4(:,ii)./t2(:,ii)*(tes(ii))^2;
    gradsum = gradsum - (1-2*prc)* t4(:,ii)*(tes(ii));%
end

hesssum(isnan(hesssum)) = 0;
% gradsum(isnan(gradsum)) = 0;
dc = -(A'*bsxfun(@times,hesssum,A))\(A'*gradsum);

% check that dc is reliable
dc(isnan(dc)) = 0;

% Update
c = c + lambda*dc;
