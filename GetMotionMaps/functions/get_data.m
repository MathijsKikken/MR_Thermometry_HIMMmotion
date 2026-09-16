function [scan_T,tes_T,GRE_info_T,file_name] = get_data(img_name_append, info_name_append, sequence, datafile, Pt, slice, rel_echoes)
%% MR thermometry of RF heating in the human brain at 7T using the HIMM MRT method with SVD-based motion correction
% get_data: function to load the dynamic data of a single subject
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   img_name_append:    specifies which image file needs to be loaded: 
%                       e.g., '' for raw data and '_registered_slice0' for 
%                       rigidly registered data
%   info_name_append:   specifies which file needs to be loaded
%   sequence:           specified whether we are looking at a baseline or
%                       heating scan
%   Pt:                 integer indicating which subject is loaded
%   slice:              usually multi-slice 2D data is acquired, this
%                       integer specified which slice is evaluated
%   rel_echoes:         list of echoes that are used for evaluations
%                       (usually the first echo is discarded)
%   rel_drift_echoes:   list of echoes that are used for evaluations of a
%                       second dataset
%                       (usually the first echo is discarded)
%
% Output: 
%   scan_T:             complex dynamic imaging data
%   tes_T:              echo times in seconds
%   GRE_info_T:         info file for the data that was loaded
%   file_name:          string representing the scan that was loaded


% Specify whether the baseline or the heated file will be loaded
if strcmp(sequence,'baseline')    % Baseline
    seq = 'MRTbase';
elseif strcmp(sequence,'heating') % Heating
    seq = 'MRTheat';
else
    disp('The parameter sequence was incorrectly specified')
end

% Find the muti-echo GRE MRT scan (and its info file) via .mat file
file_name = [datafile 'M',num2str(Pt,'%04u'),'\' seq '_mGRE_img' img_name_append '.mat'];
path_file_scan = file_name;
path_file_info = [datafile 'M',num2str(Pt,'%04u'),'\' seq '_mGRE_info' info_name_append '.mat'];

% Load data and provide a general name as it can be used in the remainder of the code
S_img = load(path_file_scan); N_img = fieldnames(S_img); eval(['GRE_scan = S_img.' N_img{1} ';']);
S_info = load(path_file_info); N_info = fieldnames(S_info); eval(['GRE_info = S_info.' N_info{1} ';']);

% The raw data has 7 dimensions, where the last 2 are GRE data
% Cast to complex and reshape the image such that there are 4 dimensions: [ nx | ny | nTE | dynamics]
if length(size(GRE_scan)) == 4
    scan = GRE_scan;
else
    scan = squeeze(GRE_scan(:,:,slice,:,:,1,1).*exp(1i*GRE_scan(:,:,slice,:,:,1,2)));
end

% Remove certain echoes: do this separately for the scan used for
% temperature measurements and the scan used for drift measurements
GRE_info_T = GRE_info;
scan_T = scan(:,:,rel_echoes,:);
GRE_info_T.imgdef.echo_time.uniq = GRE_info.imgdef.echo_time.uniq(rel_echoes);

% Acquire echo times
tes_T = GRE_info_T.imgdef.echo_time.uniq.*1e-3; 

end