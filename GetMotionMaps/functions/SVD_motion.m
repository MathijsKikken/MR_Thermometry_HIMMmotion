function [U,S,V] = SVD_motion(B0_motion_avg)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% SVD_motion: perform singular value decomposition to characterize the B0 disturbances
% in its principal components
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   B0_motion_avg:      motion fields orginating from all possible
%                       combinations of 2 dynamics within a dynamic scan
%
% Output: 
%   U:                  left singular vectors
%   S:                  singular values
%   V:                  right singular vectors


% Put all images into vectors and combine them into an array
A = reshape(B0_motion_avg,[size(B0_motion_avg,1)*size(B0_motion_avg,2),size(B0_motion_avg,3)]);

% Perform SVD
[U,S,V] = svd(A);

end