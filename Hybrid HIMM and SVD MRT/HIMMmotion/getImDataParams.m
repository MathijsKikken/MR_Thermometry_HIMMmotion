function [imDataParams] = getImDataParams(GRE_info,imgslib,file_name,frequency_system,b0,prc)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% getImDataParams: Define the image parameters and store them in a struct
% The struct is passed onto the ISMRM toolbox for WF reconstruction

if size(imgslib,4) == 1
    imDataParams.images = double(permute(imgslib,[1,2,5,4,3])); % An additional dimension is added
else
    imDataParams.images = double(permute(imgslib,[1,2,6,5,3,4])); % An additional dimension is added
end
imDataParams.filename = file_name;
imDataParams.TE = GRE_info.imgdef.echo_time.uniq(2:end).*1e-3;
imDataParams.nEchoes = length(imDataParams.TE);
imDataParams.nAcquisitions = size(imDataParams.images,5);
imDataParams.centerFreq_Hz = frequency_system;
imDataParams.FieldStrength = b0;
% imDataParams.B0dir = [0 0 1];   % subject to change
imDataParams.PrecessionIsClockwise = prc;
imDataParams.fatShift_px = str2double(GRE_info.pardef.Water_Fat_shift_pixels);
if isfield(GRE_info.imgdef,'slice_thickness_in_mm')
    imDataParams.voxelSize_mm = [GRE_info.imgdef.pixel_spacing_x_y.uniq(1) GRE_info.imgdef.pixel_spacing_x_y.uniq(2) GRE_info.imgdef.slice_thickness_in_mm.uniq];
else
   imDataParams.voxelSize_mm = [GRE_info.imgdef.pixel_spacing_x_y.uniq(1) GRE_info.imgdef.pixel_spacing_x_y.uniq(2) GRE_info.imgdef.slice_thickness.uniq];
end
imDataParams.isoCenter_REC = [size(imDataParams.images,1)/2,size(imDataParams.images,2)/2,size(imDataParams.images,3)/2];
imDataParams.readoutDimension_REC = 1;
imDataParams.readoutPolarity = [0 0];
% imDataParams.affMat_ij2xyz = [0 -0.9821 0 94.7994;0 0 -0.9000 36.4476;-0.9821 0 0 110.4910;0 0 0 1.0000];
% imDataParams.affMat_ij2raf = [0 0 0.9 -55.3524;0.9821 0 0 -94.7994;0 -0.9821 0 252.3425;0 0 0 1.0000];
% imDataParams.affMat_REC2xyz = [0 -0.9821 0 94.7994;0 0 -0.9000 36.4476;-0.9821 0 0 110.4910;0 0 0 1.0000];
% imDataParams.affMat_REC2raf = [0 0 -0.9000 36.4476;0 0.9821 0 -94.7994;-0.9821 0 0 49.6750;0 0 0 1.0000];
if GRE_info.imgdef.slice_orientation_TRA_SAG_COR.uniq == 1
    imDataParams.orientation = 'TRA';
elseif GRE_info.imgdef.slice_orientation_TRA_SAG_COR.uniq == 2
    imDataParams.orientation = 'SAG';
else
    imDataParams.orientation = 'COR';
end

end

