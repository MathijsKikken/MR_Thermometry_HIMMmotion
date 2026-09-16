function [imgs,wimgs,fimgs] = calcmeimgs(W,F,tes,b0,dw0,r2star,dww,prc,fatmodel)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% calcmeimgs: function to re-calculate/synthesize multi-echo images from provided Water
% and Fat images, given a fat spectrum model
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   W:        water image
%   F:        fat image
%   tes:      vector of tes in sec
%   b0:       field strength in Tesla
%   dw0:      off-resonance/rate map in rad/sec
%   r2star:   R2* map in /sec
%   dww:      water frequency shift due to chemical shift or temperature change in rad/sec
%   prc:      precession direction (0 or 1 for clockwise/counter clockwise)
%   fatmodel: optional fat peak model specification
%
% Outputs:
%   imgs:     complex images ([nRows nCols nEchoes]) (main output)
%   wimgs:    Water image at each echo
%   fimgs:    Fat image at each echo


if ~exist('fatmodel','var')
    fatppms = [-3.8 -3.4 -2.6 -1.94 -0.39 0.6];
    fatalphas = [0.087 0.693 0.128 0.004 0.039 0.048];
else
    fatppms = fatmodel.frequency;
    fatalphas = fatmodel.relAmps;
end

fatfrqs = b0*42.577481*fatppms; % fat offsets, Hz

wimgs = zeros([size(W) length(tes)]);
fimgs = wimgs;
for ii = 1:length(tes)
    % water contribution
    wimgs(:,:,ii) = W.*exp(-(1-2*prc)*1i*dww*tes(ii));%
    % off-resonance/rate map
    wimgs(:,:,ii) = wimgs(:,:,ii).*exp((-(1-2*prc)*1i*dw0-r2star)*tes(ii)); %
    % fat contributions
    for jj = 1:length(fatfrqs)
        fimgs(:,:,ii) = fimgs(:,:,ii) + F.*fatalphas(jj)*exp(-(1-2*prc)*1i*2*pi*fatfrqs(jj)*tes(ii));% For some reason, making this prc -(1-2*prc) and the others (1-2*prc) works for simulation but not dog
    end
    % off-resonance/rate map
    fimgs(:,:,ii) = fimgs(:,:,ii).*exp((-(1-2*prc)*1i*dw0-r2star)*tes(ii)); %
    
end
imgs = wimgs + fimgs;