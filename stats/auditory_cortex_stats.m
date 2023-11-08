%%%
%   NeuroTin / 20231101 / auditory_cortex_stats.m
%
%   Down-reg. eval. of auditory cortex ROIs.
%%%
%
%   Nicolas Gninenko / nicolas.gninenko@gmail.com
%
%%%


% Load the spmT map
tmp_auditory = load_untouch_nii('spmT_0001_maskedwith_AudROI_LR_21subj_Imask_MNI.nii');
tmp_auditory = double(tmp_auditory.img);


% x = 46 is separating well the hemispheres (left/right)
% but pay attention that it is in RADIOLOGIC CONVENTION,
% i.e., right is left and vice versa
tmp_auditory_L = tmp_auditory(46:end,:,:);
tmp_auditory_R = tmp_auditory(1:45,:,:); % visual sanity check -> OK

size_ROI_L = length(find(tmp_auditory_L));
size_ROI_R = length(find(tmp_auditory_R));
max_val_L = max(max(max(tmp_auditory_L))); max_val_R = max(max(max(tmp_auditory_R)));
min_val_L = min(min(min(tmp_auditory_L))); min_val_R = min(min(min(tmp_auditory_R)));


% Show mean ± [SD] for both sides
fprintf(['\nAverage down-regulation (left: mean ± [SD]): ' num2str(mean(nonzeros(tmp_auditory_L(:)))) ...
    ' ± ' num2str(std(nonzeros(tmp_auditory_L(:)))) ' (' num2str(length(nonzeros(tmp_auditory_L))) ' voxels)\n']);
fprintf(['Average down-regulation (right: mean ± [SD]): ' num2str(mean(nonzeros(tmp_auditory_R(:)))) ...
    ' ± ' num2str(std(nonzeros(tmp_auditory_R(:)))) ' (' num2str(length(nonzeros(tmp_auditory_R))) ' voxels)\n']);



