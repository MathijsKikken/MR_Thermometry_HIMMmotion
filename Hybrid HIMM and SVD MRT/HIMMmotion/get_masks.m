function [mask_F,body_mask,brain_mask] = get_masks(scan, body_factor, body_mask_mGRE)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% get_masks: function to create the masks used in the HIMMmotion algorithm
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   scan:               gradient echo scan of serveral echo times
%   body_factor:        factor used to define the body mask
%                       
% Output: 
%   mask_F:             fat mask
%   body_mask:          mask of the regions that have tissues
%   brain_mask:         brain mask
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify region that can be seen as being human body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Thresholding based on magnitude image of the first GRE (and closing)
    body_mask = zeros([size(scan,1) size(scan,2)]);
    body_mask(abs(scan(:,:,1)) > body_factor*mean(abs(scan(:,:,1)),'all')) = 1;
    body_mask = imclose(body_mask,strel("disk",10));
    body_mask = body_mask.*body_mask_mGRE;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the fat map by evaluating the gradient image of the first echo magnitude image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Use otsu thresholding on the gradient of the first echo magnitude image
    img_temp = abs(scan(:,:,1));
    img_temp_body = img_temp(body_mask == 1);
    mean_img = mean(img_temp_body);
    img_temp(img_temp > mean_img) = mean_img;

    [Gmag,~] = imgradient(img_temp);
    Gmag_norm = mat2gray(Gmag);
    level_gradient = graythresh(Gmag_norm(body_mask == 1));
    mask_F = zeros(size(Gmag));
    mask_F(Gmag_norm > level_gradient) = 1;
    mask_F = imclose(mask_F,strel('disk',15));
    mask_F = mask_F.*body_mask;    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the brain mask from the body mask and fat mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    brain_mask = body_mask-mask_F;
    brain_mask(brain_mask < 0) = 0;
    
end