%%%
%   NeuroTin / 20201007 / ntin_preproc_NF_new.m
%
%   Script to batch pre-process any NF data from a specified NeuroTin
%   subject. Steps are described below in more details.
%   Adapted for running on srv3 and/or srv4 (MIP:Lab).
%
%
%
%   Usage : ntin_preproc_NF_new(subj_ID, urflag, FD_th, sflag);
%
%   - where 'subj_ID' is for example '002GNV' (subject ID)
%   - where 'urflag' is either 0 or 1, 0 for gre_field_unwarping
%                                      1 for normal realign & reslice
%   - where 'FD_th' is a numeric value for FD thresholding
%   - where 'sflag' is either 0 or 1, 0 for not smoothing anything
%                                     1 for smoothing Reg* volumes, etc.
%
%
%
%   Pre-processing steps are (steps are automatically skipped if needed) :
%
%
%   0) Create ../Processed/000GNV/[NF_dir] and its /_misc/ directory;
%       then import all NF raw DICOMs (to .nii) into the above folder;          =======> f*.nii
%   [0.1) The ../Data/[all_00GNV_Vall]/_posthoc/gre_field_nii dirs are
%       already created from the ntin_preproc_rsfMRI_new script; that
%       should be run prior to neurofeedback pre-processing! So that the
%       phase & magnitude maps are already available for NF data if needed]
%
%   1) Slice timing correction, with following parameters:                      =======> af*.nii
%       - TR 1.5 [s]
%       - TA = 1.5 - (1.5/64) [s]
%       - nslices = 64
%       - slice order = [2:2:64 1:2:63]
%       - ref slice = 2 (i.e. first slice, by default)
%
%   2) Realign & unwarp (uf*) OR/AND realign & reslice (rf*) all inputs
%       according to the subject-specific suggested pre-processing flag;        =======> mean*.nii (meanaf* or meanuaf*)
%   3) Coregister all realigned/unwarped volumes to T1 V01 anatomical;          =======> uaf*.nii or raf*.nii, as well as rp_*.txt motion params file
%   4) Create intersect_mask_NF_000GNV_rf/uf.nii from all vols;                 =======> intersect_mask_NF_000GNV_raf/uaf.nii
%   4.1) Bet the intersect_mask_NF_000GNV_rf/uf.nii (slightly) in order
%       to obtain a cleaner masked volume as well as its mask;                  =======> intersect_mask_NF_000GNV_raf/uaf_betted(_mask).nii
%   4.2) Bet the mean*.nii template from stage 1);                              =======> mmean*.nii & mmean*_mask.nii
%
%   5) Run PhysIO TAPAS in non-verbose mode in order to create the
%       regressors of interest related to motion (12, including
%       derivatives), physiological noise (16, from cardiac & breathing
%       belt), and FDs as specified from input (FD_th).
%   - Note that FDs are not regressed out but only kept for information;        =======> multiple_regressors.txt (for each visit)
%
%   6) Run T1_V01 segmentation & CSF/WM masks remapping for further
%       regression in the following step (i.e. creates the T1_V01 folder
%       and runs the T1 V01 segmentation if not already done; then maps
%       back the CSF & WM ROIs to functional space according to the
%       template obtained at step 1));                                          =======> irwc2WhiteMask_09_108x108x64_rsfMRI.nii
%                                                                                        irwc3CsfMask_07_108x108x64_rsfMRI.nii
%
%   7) Proper regression of nuisance variables as mentioned above;              =======> Reg2_Physio_*.nii and/or Reg2_MotionOnly_*.nii
%
%   8) ROIs unwarping/realignment for further analysis;                         =======> remapped_ROIs folder in NF_dir (preprocessed directory)
%
%   9) Smoothing (optional) step for Reg2_* files above;                        =======> s6Reg2_*.nii
%   9.1) Intersect mask creation for smoothed volumes;                          =======> intersect_mask_NF_000GNV_s6Reg2.nii
%
%   [[10) (Optional) normalization to MNI step for Reg* files aboves;
%   =======> wReg*.nii]]
%%%
%
%   Nicolas Gninenko / nicolas.gninenko@gmail.com
%
%%%

function ntin_preproc_NF_new(subj_ID)




%   Default [and some subject-specific] parameters

%subj_ID = '002GNV';

aflag   = 1; % carry out slice timing correction
which_visit_as_r_u_template = 'V01';
physiological_folder_name = '_physiological_cardiac';
reg_vols_prefix = 'Reg2_'; % '2' is for preprocessing pipeline consistency with NF only
FD_th   = 0.5;
urflag  = 0;
fsl_bet_default_th_avgvol = '0.36'; % for average volume (mean*.nii)
fsl_bet_default_th_intvol = '0.2'; % for intersect mask (intersect_mask*.nii)
sflag   = 1; % also create smoothed volumes (from the Reg2_*) by default
smoothing_factor = 6;
nflag   = 1; % normalize to MNI space using T1_V01 segmentation's warping?
which_regression_approach = 2; % use suggested modifications by E. Amico regarding motion (no FDs!), physiological noise, and CSF/WM regression
switch subj_ID
    case '002GNV'
        %urflag = 1; % unwarping used // not used for final group analysis
        fsl_bet_default_th_avgvol = '0.47'; % try a higher threshold for 002GNV because of strong frontal artefacts
    case '003GNV'
        %which_visit_as_r_u_template = 'V04'; % not used
        %urflag = 1; % unwarping used // same as above
        % 003GNV: V11 Run4 is corrupted -> has been moved to _posthoc/ in the corresponding raw (Data/) session folder
        fsl_bet_default_th_avgvol = '0.25'; %
    case '004GNV'
        %
        fsl_bet_default_th_avgvol = '0.25';
    case '005GNV'
        which_visit_as_r_u_template = 'V03';
        %urflag = 1; % trying without
    case '008GNV'
        FD_th = 0.7;
        which_visit_as_r_u_template = 'V09';
        %urflag = 1; % trying without
    case '012GNV' % done
        urflag = 0;
        which_visit_as_r_u_template = 'V07';
        FD_th = 0.7;
    case '013GNV'
        % 013GNV: V13 gre_f map is very likely corrupted -> do not use
    case '016GNV' % done
        %urflag = 0;
        which_visit_as_r_u_template = 'V03';
    case '020GNV'
        fsl_bet_default_th_avgvol = '0.5';
    case '021GNV'
        urflag = 0;
        FD_th = 0.6;
        which_visit_as_r_u_template = 'V04';
    case '023GNV'
        urflag = 0;
        FD_th = 0.6;
    case '025GNV'
        urflag = 0;
    case '026GNV'
        which_visit_as_r_u_template = 'V02';
        fsl_bet_default_th_avgvol = '0.38';
    case '034GNV'
        urflag = 0;
        %FD_th = 0.5;
        fsl_bet_default_th_avgvol = '0.42';
    case '037GNV'
        urflag = 0;
        which_visit_as_r_u_template = 'V05';
        fsl_bet_default_th_avgvol = '0.4';
    case '041GNV'
        urflag = 0;
        fsl_bet_default_th_avgvol = '0.4';
    case '043GNV'
        urflag = 0;
        which_visit_as_r_u_template = 'V09';
    case '048GNV'
        urflag = 0;
        which_visit_as_r_u_template = 'V10';
    case '049GNV'
        fsl_bet_default_th_avgvol = '0.3';
    case '054GNV'
        urflag = 0;
    case '056GNV'
        urflag = 0;
        fsl_bet_default_th_avgvol = '0.3';
        which_visit_as_r_u_template = 'V04';
    otherwise
        error(['Unknown subject specified (' subj_ID ')... Please double-check your input!']);
end
NF_dir = 'NF_new';

if urflag==1
    r_u_prefix = 'u'; % uf*.nii outputs
elseif urflag==0
    r_u_prefix = 'r'; % rf*.nii outputs
end
if aflag==1
    aflag_prefix = 'a';
elseif aflag==0
    aflag_prefix = '';
end




%   General paths & vars

if ismac && strcmp(getenv('USER'),'nicogn')
    basepath_data = [filesep 'Volumes' filesep 'NeuroTin2' filesep 'Data' filesep];
    procpath = [filesep 'Volumes' filesep 'NeuroTin2' filesep 'Processed' filesep];
    spm_dir = [filesep 'Users' filesep 'nicogn' filesep 'Documents' filesep 'MATLAB' filesep 'spm12' filesep];
    % T1 V01 part (tmp) below
    addpath('/Volumes/NeuroTin2/MATLAB/NIfTI_20140122/');
    addpath(genpath('/Volumes/NeuroTin2/MATLAB/static_FC_pipeline/DPARSF_V2.3_130615/'));
    addpath(genpath('/Volumes/NeuroTin2/MATLAB/static_FC_pipeline/depMat/'));
    addpath('/Volumes/NeuroTin2/MATLAB/static_FC_pipeline/depMat_addon/');
    addpath('/Volumes/NeuroTin2/MATLAB/static_FC_pipeline/preproc_SPM8_original_qnap/');
elseif isunix && (strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv4') || ...
        strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv3') || ...
        strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv2') || ...
        strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv1'))
    %
    % this part of the code (paths on internal servers) is not published for privacy reasons
    % [adapt paths above to use the code]
    %
    %procpath = [];
    tmp_template_EPI_file_path = [procpath 'Templates' filesep 'example_template_EPI_Aud_NF_Run_seq.dcm'];
end
addpath(spm_dir);
addpath([spm_dir 'matlabbatch' filesep]);
%spm_dflt_gre_f_templ_mask_img = [spm_dir 'toolbox' filesep 'FieldMap' filesep 'T1.nii'];
templates_wmcsf_path = [procpath 'Templates' filesep 'GM_CSF' filesep];
fsl_bin_path = [filesep 'usr' filesep 'local' filesep 'fsl' filesep 'bin' filesep]; % common to local & srv settings
setenv('FSLDIR',[filesep 'usr' filesep 'local' filesep 'fsl']);
setenv('FSLOUTPUTTYPE','NIFTI');
func_vols_size = [108,108,64];
CSF_template_name = 'CsfMask_07_121x145x121.nii';
WM_template_name = 'WhiteMask_09_121x145x121.nii';
CSF_template_name_func = ['CsfMask_07_' num2str(func_vols_size(1)) 'x' num2str(func_vols_size(2)) 'x' num2str(func_vols_size(3)) '.nii']; % (corresp to c3)
WM_template_name_func = ['WhiteMask_09_' num2str(func_vols_size(1)) 'x' num2str(func_vols_size(2)) 'x' num2str(func_vols_size(3)) '.nii']; % (corresp to c2)

tmp_ini_pwd = pwd;




%   Parallel computing params for regression

if isunix && strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv4')
    num_workers_parallel = 72; % srv 4
elseif isunix && strcmp(char(java.net.InetAddress.getLocalHost.getHostName),'miplabsrv3')
    num_workers_parallel = 56; % srv 3
elseif isunix
    num_workers_parallel = 24; % srv 2 or 1 or other
elseif ismac
    num_workers_parallel = 0; % don't use locally
end




%   0) Create ../Processed/000GNV/[NF_dir] and its /_misc/ directory;
%   	then import all rsfMRI raw DICOMs (to .nii) into the above folder;

all_data_dirs = dir([basepath_data '*' subj_ID '_V*']);
if ~isequal(length(all_data_dirs),15), warning(['Unable to find all the visits for participant ' subj_ID '! Check manually...']); end
if ~exist([procpath subj_ID filesep NF_dir filesep],'dir'), mkdir([procpath subj_ID filesep NF_dir filesep]); end
if ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep],'dir'), mkdir([procpath subj_ID filesep NF_dir filesep '_misc' filesep]); end
c = 1; clear matlabbatch;
for v = 1:length(all_data_dirs)
    %if v<10, v_str = ['V0' num2str(v)]; else, v_str = ['V' num2str(v)]; end
    v_str = all_data_dirs(v).name((regexp(all_data_dirs(v).name,'GNV_V')+4):(regexp(all_data_dirs(v).name,'GNV_V')+6)); % find v_str from file name instead of assuming from v value
    if ~exist([procpath subj_ID filesep NF_dir filesep v_str filesep],'dir')
        mkdir([procpath subj_ID filesep NF_dir filesep v_str filesep]);
    end
    all_runs_per_session = dir([basepath_data all_data_dirs(v).name filesep '*_Aud_NF_Run_*']); % find all the NF runs
    for h = 1:length(all_runs_per_session)
        clear tmp_runfoldername;
        tmp_run_foldername = ['Run' num2str(all_runs_per_session(h).name(regexp(all_runs_per_session(h).name,'Run_')+4))]; % Run0 by default, as output folder name
        if ~isempty(regexp(all_runs_per_session(h).name,'Transfer','once')), tmp_run_foldername = [tmp_run_foldername '_Transfer']; %#ok
        elseif ~isempty(regexp(all_runs_per_session(h).name,'NoFeedback','once')), tmp_run_foldername = [tmp_run_foldername '_NoFeedback']; %#ok
        end
        if ~exist([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_run_foldername filesep],'dir')
            mkdir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_run_foldername filesep]);
        elseif isequal(length(dir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_run_foldername filesep 'f*.nii'])),270)
            fprintf(['\nSkipping NF DICOM import for ' subj_ID ' ' v_str ' ' tmp_run_foldername ' (files already exist)...\n']);
            continue;
        end
        all_raw_data_per_run = dir([basepath_data all_data_dirs(v).name filesep all_runs_per_session(h).name filesep 'raw' filesep '*.dcm']);
        if ~isequal(length(all_raw_data_per_run),270)
            error(['Couldn''t find all the NF volumes (only ' num2str(length(all_raw_data_per_run)) '/270) for ' all_data_dirs(v).name ', ' all_runs_per_session(h).name '... Please check manually.']);
            %continue; % skipping
        end
        matlabbatch{c}.spm.util.import.dicom.data = cell(length(all_raw_data_per_run),1); %#ok
        for hh = 1:length(all_raw_data_per_run)
            matlabbatch{c}.spm.util.import.dicom.data{hh,1} = ...
                [basepath_data all_data_dirs(v).name filesep all_runs_per_session(h).name filesep 'raw' filesep all_raw_data_per_run(hh).name];
        end
        matlabbatch{c}.spm.util.import.dicom.outdir{1} = [procpath subj_ID filesep NF_dir filesep v_str filesep tmp_run_foldername filesep]; %#ok
        matlabbatch{c}.spm.util.import.dicom.root = 'flat'; %#ok
        matlabbatch{c}.spm.util.import.dicom.protfilter = '.*'; %#ok
        matlabbatch{c}.spm.util.import.dicom.convopts.format = 'nii'; %#ok
        matlabbatch{c}.spm.util.import.dicom.convopts.icedims = 0; %#ok
        if ~exist(matlabbatch{c}.spm.util.import.dicom.outdir{1},'dir'), mkdir(matlabbatch{c}.spm.util.import.dicom.outdir{1}); end
        c = c+1;
    end
end
if c>1
    fprintf(['\nImporting all (or remaining) NF DICOMs for ' subj_ID '...\n']); tic;
    %spm('defaults','fmri');
    %spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
    clear matlabbatch c tmp_run_foldername all_raw_data_per_run all_runs_per_session;
end




%   1) Slice timing correction;

clear matlabbatch;
c = 1;
if aflag
    if ~isempty(dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep aflag_prefix 'f*.nii']))
        fprintf(['\nSkipping slice timing correction step for ' subj_ID ' (' aflag_prefix 'f*.nii volumes already exist in the V01 Run1 folder)...\n']);
    else
        for v = 1:length(all_data_dirs)
            v_str = all_data_dirs(v).name((regexp(all_data_dirs(v).name,'GNV_V')+4):(regexp(all_data_dirs(v).name,'GNV_V')+6)); % find v_str from file name instead of assuming from v value
            all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep v_str filesep 'Run*']);
            for r = 1:length(all_NF_runs_per_visit)
                tmp_all_scans = dir([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_visit(r).name filesep 'f*.nii']);
                if length(tmp_all_scans)~=270, error(['Invalid number of f*.nii inputs found for ' subj_ID ' ' v_str ' ' all_NF_runs_per_visit(r).name '. Please check manually or re-run DICOM import.']); end
                matlabbatch{1}.spm.temporal.st.scans{1,c} = cell(length(tmp_all_scans),1);
                for rr = 1:length(tmp_all_scans)
                    matlabbatch{1}.spm.temporal.st.scans{1,c}{rr,1} = [procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_visit(r).name filesep tmp_all_scans(rr).name ',1'];
                end
                c = c+1;
            end
        end
        tmp_so = dicominfo(tmp_template_EPI_file_path);
        matlabbatch{1}.spm.temporal.st.nslices = 64;
        matlabbatch{1}.spm.temporal.st.tr = 1.5;
        matlabbatch{1}.spm.temporal.st.ta = 1.5-(1.5/64); %1.4765625;
        matlabbatch{1}.spm.temporal.st.so = transpose(typecast(tmp_so.Private_0019_1029,'double'));
        matlabbatch{1}.spm.temporal.st.refslice = 0; %2; %matlabbatch{1}.spm.temporal.st.so(1);
        matlabbatch{1}.spm.temporal.st.prefix = aflag_prefix; % 'a' usually
        fprintf(['\nRunning slice timing correction (NF) for ' subj_ID '...\n']); tic;
        %spm('defaults','fmri');
        %spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
        clear matlabbatch c tmp_all_scans all_NF_runs_per_visit tmp_so;
    end
else
    fprintf(['\nSkipping slice timing correction for ' subj_ID ' (because aflag = 0).\n']);
end




%   2) Realign & unwarp (uf*) OR/AND realign & reslice (rf*) all inputs
%       according to the subject-specific suggested pre-processing flag;

tmp_skip_coreg = 0;
clear matlabbatch;
c = 1;
if urflag==1
    if isequal(length(dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep r_u_prefix aflag_prefix '*.nii'])),270)
        fprintf(['\nSkipping realign & unwarp step for ' subj_ID ' (' r_u_prefix aflag_prefix '*.nii volumes already exist in the V01 Run1 folder)...\n']);
        tmp_skip_coreg = 1;
    else
        cd([procpath subj_ID filesep NF_dir filesep '_misc' filesep]);
        tmp_nb_of_runs_per_visit = zeros(length(all_data_dirs),1);
        for v = 1:length(all_data_dirs)
            v_str = all_data_dirs(v).name((regexp(all_data_dirs(v).name,'GNV_V')+4):(regexp(all_data_dirs(v).name,'GNV_V')+6)); % find v_str from file name instead of assuming from v value
            tmp_all_runs = dir([procpath subj_ID filesep NF_dir filesep v_str filesep 'Run*']);
            tmp_nb_of_runs_per_visit(v,1) = length(tmp_all_runs);
            for r = 1:tmp_nb_of_runs_per_visit(v,1)
                %   if ~isempty(dir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep 'u*.nii'])), continue; end
                tmp_all_scans = dir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep aflag_prefix 'f*.nii']);
                if ~isequal(length(tmp_all_scans),270)
                    error(['Unable to find the correct amount of imported DICOMs (or slice-timing corrected volumes; found ' num2str(length(tmp_all_scans)) '/270) for ' subj_ID ' ' v_str ' ' tmp_all_runs(r).name '... Please check manually.']);
                end
                for z = 1:length(tmp_all_scans)
                    matlabbatch{1}.spm.spatial.realignunwarp.data(c).scans{z,1} = ...
                        [procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep tmp_all_scans(z).name ',1'];
                end
                tmp_vdm_path = dir([basepath_data all_data_dirs(v).name filesep '_posthoc' filesep 'gre_field_nii' filesep 'vdm5*.nii']);
                if isempty(tmp_vdm_path), error(['Unable to find VDM map for session ' all_data_dirs(v).name ' for ' subj_ID '! Please check manually.']); end
                matlabbatch{1}.spm.spatial.realignunwarp.data(c).pmscan{1} = ...
                    [basepath_data all_data_dirs(v).name filesep '_posthoc' filesep 'gre_field_nii' filesep tmp_vdm_path(1).name ',1'];
                c = c+1;
            end
        end
        % for some participants, we use a different [better] template
        % for realignment
        if ~strcmp(which_visit_as_r_u_template,'V01')
            if strcmp(subj_ID,'005GNV') && strcmp(which_visit_as_r_u_template,'V03')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(13);
                matlabbatch{1}.spm.spatial.realignunwarp.data(13) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'008GNV') && strcmp(which_visit_as_r_u_template,'V09')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(57);
                matlabbatch{1}.spm.spatial.realignunwarp.data(57) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'012GNV') && strcmp(which_visit_as_r_u_template,'V07')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(44);
                matlabbatch{1}.spm.spatial.realignunwarp.data(44) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'016GNV') && strcmp(which_visit_as_r_u_template,'V03')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(15);
                matlabbatch{1}.spm.spatial.realignunwarp.data(15) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'021GNV') && strcmp(which_visit_as_r_u_template,'V04')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(22);
                matlabbatch{1}.spm.spatial.realignunwarp.data(22) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'026GNV') && strcmp(which_visit_as_r_u_template,'V02')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(7);
                matlabbatch{1}.spm.spatial.realignunwarp.data(7) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'037GNV') && strcmp(which_visit_as_r_u_template,'V05')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(25);
                matlabbatch{1}.spm.spatial.realignunwarp.data(25) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'043GNV') && strcmp(which_visit_as_r_u_template,'V09')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(49);
                matlabbatch{1}.spm.spatial.realignunwarp.data(49) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'048GNV') && strcmp(which_visit_as_r_u_template,'V10')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(55);
                matlabbatch{1}.spm.spatial.realignunwarp.data(55) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'056GNV') && strcmp(which_visit_as_r_u_template,'V04')
                tmp_ttt = matlabbatch{1}.spm.spatial.realignunwarp.data(1);
                matlabbatch{1}.spm.spatial.realignunwarp.data(1)  = matlabbatch{1}.spm.spatial.realignunwarp.data(19);
                matlabbatch{1}.spm.spatial.realignunwarp.data(19) = tmp_ttt;
                clear tmp_ttt;
            end
        end
        %
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0; % do not do 2nd pass, it is ok
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 0; % do not mask as well, it is ok
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = r_u_prefix; % 'u' or 'r' for NF new pre-proc pipeline
        %
        fprintf(['\nRunning realign & unwarp (NF, ' which_visit_as_r_u_template ') for ' subj_ID '...\n']); tic;
        spm_jobman('run',matlabbatch);
        fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
        cd(tmp_ini_pwd);
        clear matlabbatch tmp_vdm_path tmp_all_runs tmp_all_scans tmp_nb_of_runs_per_visit; % do not clear tmp_skip_coreg var
    end
elseif urflag==0
    if isequal(length(dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep r_u_prefix aflag_prefix '*.nii'])),270)
        fprintf(['\nSkipping realign & reslice step for ' subj_ID ' (' r_u_prefix aflag_prefix '*.nii volumes already exist in the V01 Run1 folder)...\n']);
        tmp_skip_coreg = 1;
    else
        cd([procpath subj_ID filesep NF_dir filesep '_misc' filesep]);
        for v = 1:length(all_data_dirs)
            %if v<10, v_str = ['V0' num2str(v)]; else, v_str = ['V' num2str(v)]; end
            v_str = all_data_dirs(v).name((regexp(all_data_dirs(v).name,'GNV_V')+4):(regexp(all_data_dirs(v).name,'GNV_V')+6)); % find v_str from file name instead of assuming from v value
            tmp_all_runs = dir([procpath subj_ID filesep NF_dir filesep v_str filesep 'Run*']);
            for r = 1:length(tmp_all_runs)
                %   if ~isempty(dir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep 'r*.nii'])), continue; end
                tmp_all_scans = dir([procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep aflag_prefix 'f*.nii']);
                if ~isequal(length(tmp_all_scans),270)
                    error(['Unable to find the correct amount of imported DICOMs (or slice-timing corrected volumes; found ' num2str(length(tmp_all_scans)) '/270) for ' subj_ID ' ' v_str ' ' tmp_all_runs(r).name '... Please check manually.']);
                end
                matlabbatch{1}.spm.spatial.realign.estwrite.data{1,c} = cell(270,1);
                for z = 1:length(tmp_all_scans)
                    matlabbatch{1}.spm.spatial.realign.estwrite.data{1,c}{z,1} = ...
                        [procpath subj_ID filesep NF_dir filesep v_str filesep tmp_all_runs(r).name filesep tmp_all_scans(z).name ',1'];
                end
                c = c+1;
            end
        end
        % same as previously
        if ~strcmp(which_visit_as_r_u_template,'V01')
            if strcmp(subj_ID,'005GNV') && strcmp(which_visit_as_r_u_template,'V03')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(13);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(13) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'008GNV') && strcmp(which_visit_as_r_u_template,'V09')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(57);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(57) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'012GNV') && strcmp(which_visit_as_r_u_template,'V07')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(44);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(44) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'016GNV') && strcmp(which_visit_as_r_u_template,'V03')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(15);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(15) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'021GNV') && strcmp(which_visit_as_r_u_template,'V04')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(22);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(22) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'026GNV') && strcmp(which_visit_as_r_u_template,'V02')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(7);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(7) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'037GNV') && strcmp(which_visit_as_r_u_template,'V05')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(25);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(25) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'043GNV') && strcmp(which_visit_as_r_u_template,'V09')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(49);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(49) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'048GNV') && strcmp(which_visit_as_r_u_template,'V10')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(55);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(55) = tmp_ttt;
                clear tmp_ttt;
            elseif strcmp(subj_ID,'056GNV') && strcmp(which_visit_as_r_u_template,'V04')
                tmp_ttt = matlabbatch{1}.spm.spatial.realign.estwrite.data(1);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(1) = matlabbatch{1}.spm.spatial.realign.estwrite.data(19);
                matlabbatch{1}.spm.spatial.realign.estwrite.data(19) = tmp_ttt;
                clear tmp_ttt;
            end
        end
        %
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = r_u_prefix;
        %
        fprintf(['\nRunning realign & reslice (NF, ' which_visit_as_r_u_template ') for ' subj_ID '...\n']); tic;
        spm_jobman('run',matlabbatch);
        fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
        cd(tmp_ini_pwd);
        clear matlabbatch tmp_all_runs tmp_all_scans; % do not clear tmp_skip_coreg var!
    end
end




%   3) Coregister all realigned/unwarped volumes to T1 V01 anatomical;

if ~tmp_skip_coreg % for now we always coreg to T1 V01
    %clear matlabbatch;
    tmp_T1_path = dir([basepath_data all_data_dirs(1).name filesep 'T1_D2' filesep 's*.nii']);
    if isempty(tmp_T1_path), error(['Unable to find anatomical T1 volume for V01 for ' subj_ID ' for coregistration! Please double-check manually...']); end
    cd([procpath subj_ID filesep NF_dir filesep '_misc' filesep]);
    %if ~strcmp(which_visit_as_r_u_template,'V01') % [just to double-check if Run1 is named differently (e.g. Run1_Transfer) on an extremely rare occasion?]
    tmp_run_name = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*']);
    if isempty(tmp_run_name)
        error(['Unable to find the preprocessed Run1 folder for ' subj_ID ' ' which_visit_as_r_u_template '; maybe it does not exist for this visit? Please check manually...']);
    end
    tmp_run_name = tmp_run_name(1).name; % [just to double-check if Run1 is named differently (e.g. Run1_Transfer) on an extremely rare occasion?]
    %else, tmp_run_name = 'Run1'; end
    if urflag==1
        tmp_mean_source = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean' r_u_prefix aflag_prefix 'f*.nii']); % meanu(a)f*.nii
    elseif urflag==0
        tmp_mean_source = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean' aflag_prefix 'f*.nii']);
    end
    if isempty(tmp_mean_source), error(['Unable to find (' r_u_prefix aflag_prefix ')mean EPI template for ' subj_ID ' for coregistration! Please double-check manually...']); end
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1} = [basepath_data all_data_dirs(1).name filesep 'T1_D2' filesep tmp_T1_path(1).name ',1'];
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = [procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep tmp_mean_source(1).name ',1'];
    all_NF_visits = dir([procpath subj_ID filesep NF_dir filesep 'V*']);
    c = 1;
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            tmp_all_scans = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep r_u_prefix aflag_prefix 'f*.nii']);
            if length(tmp_all_scans)~=270
                error(['Incorrect number of functional (' aflag_prefix 'f*.nii) volumes found during coregistration for ' subj_ID ' (' num2str(length(tmp_all_scans)) '/270)! Please check manually.']);
            end
            for f = 1:length(tmp_all_scans)
                matlabbatch{1}.spm.spatial.coreg.estimate.other{c,1} = ...
                    [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep tmp_all_scans(f).name ',1'];
                c = c+1;
            end
        end
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    %
    fprintf(['\nRunning coregistration (NF) for ' subj_ID '...\n']); tic;
    spm_jobman('run',matlabbatch);
    fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
    cd(tmp_ini_pwd);
    clear matlabbatch tmp_all_scans tmp_mean_source tmp_T1_path tmp_run_name all_NF_runs_per_visit all_NF_visits;
else
    fprintf(['\nSkipping coregistration (since previous step already carried out - either realign or r&unwarp...) for ' subj_ID '...\n']);
end




%   4) Create intersect_mask_NF_000GNV_rf/uf.nii from all vols;

if ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f.nii'],'file')
    fprintf(['\nRunning intersect mask creation (NF, from ' r_u_prefix aflag_prefix 'f*.nii volumes) for ' subj_ID '...\n']); tic;
    all_NF_visits = dir([procpath subj_ID filesep NF_dir filesep 'V*']);
    tmp_run_name = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*']);
    if isempty(tmp_run_name)
        error(['Unable to find the preprocessed Run1 folder for ' subj_ID ' ' which_visit_as_r_u_template '; maybe it does not exist for this visit? Please check manually...']);
    end
    tmp_run_name = tmp_run_name(1).name;
    intersect_mask = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean*.nii']); % does not matter whether meanuf/meanf
    if isempty(intersect_mask)
        error(['Unable to find (' r_u_prefix aflag_prefix ')mean EPI template for ' subj_ID ' (' which_visit_as_r_u_template ' ' ...
            tmp_run_name ') for intersect mask creation... Please double-check manually...']);
    end
    intersect_mask = load_untouch_nii([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep intersect_mask(1).name]);
    intersect_mask_data = double(intersect_mask.img);
    intersect_mask_data = intersect_mask_data(:);
    intersect_mask.fileprefix = ['intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f'];
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            all_urf_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep r_u_prefix aflag_prefix 'f*.nii']);
            if ~isequal(length(all_urf_files),270)
                error(['Incorrect number of ' r_u_prefix aflag_prefix 'f*.nii volumes found in folder ' all_NF_visits(v).name ' ' ...
                    all_NF_runs_per_visit(r).name ' for ' subj_ID ' during intersect mask creation. Please check manually!']);
            end
            for f = 1:length(all_urf_files)
                tmp_vol_for_erosion = load_untouch_nii(...
                    [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_urf_files(f).name]);
                tmp_vol_for_erosion = reshape(double(tmp_vol_for_erosion.img),size(tmp_vol_for_erosion.img,1)*...
                    size(tmp_vol_for_erosion.img,2)*size(tmp_vol_for_erosion.img,3),1);
                %intersect_mask_data = max(min(),0); % max precaution?
                intersect_mask_data = min(intersect_mask_data,tmp_vol_for_erosion);
            end
        end
    end
    intersect_mask.img = int16(reshape(intersect_mask_data,size(tmp_vol_for_erosion,1),size(tmp_vol_for_erosion,2),size(tmp_vol_for_erosion,3)));
    save_untouch_nii(intersect_mask,[procpath subj_ID filesep NF_dir filesep '_misc' filesep intersect_mask.fileprefix '.nii']);
    fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
    %cd(tmp_pwd);
    clear matlabbatch all_scans_tmp tmp_mean_source tmp_T1_path;
    clear intersect_mask* tmp_vol_for_erosion all_urf_files all_NF_visits all_NF_runs_per_visit;
    fprintf(['\nDone! (' num2str(toc) ' s)\n']);
else
    fprintf(['\nSkipping intersect mask creation (already exists) for ' subj_ID '...\n']);
end




%   4.1) Bet the intersect_mask_NF_000GNV_rf/uf.nii (slightly) in order
%       to obtain a cleaner masked volume as well as its mask;

if ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f_betted.nii'],'file') ...
        || ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f_betted_mask.nii'],'file')
    fprintf(['\nBetting (f = ' fsl_bet_default_th_intvol ') intersect mask (NF) for ' subj_ID '...\n']); tic;
    [~,tmp_msg] = system([fsl_bin_path 'bet ' procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f.nii' ...
                ' ' procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' r_u_prefix aflag_prefix 'f_betted.nii' ...
                ' -f ' fsl_bet_default_th_intvol ' -m']); % default options for now
    if ~isempty(tmp_msg), disp(tmp_msg); end
    clear tmp_msg;
    fprintf(['\nDone! (' num2str(toc) ' s)\n']);
else
    fprintf(['\nSkipping intersect mask betting (already exists) for ' subj_ID '...\n']);
end




%   4.2) Bet the mean*.nii template from stage 1);

tmp_run_name = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*']);
if isempty(tmp_run_name)
    error(['Unable to find the preprocessed Run1 folder for ' subj_ID ' ' which_visit_as_r_u_template '; maybe it does not exist for this visit? Please check manually...']);
end
tmp_run_name = tmp_run_name(1).name;
if isempty(dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mmean*.nii'])) ...
        || isempty(dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mmean*_mask.nii']))
    fprintf(['\nBetting (f = ' fsl_bet_default_th_avgvol ') mean*.nii volume (NF) for ' subj_ID ' ' which_visit_as_r_u_template ' ' tmp_run_name '...\n']); tic;
    tmp_unbetted_meanurf = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean*.nii']);
    if isempty(tmp_unbetted_meanurf), error(['Unable to find mean*.nii volume (NF) for ' subj_ID ' ' which_visit_as_r_u_template ' ' tmp_run_name '... Please double-check!']); end
    tmp_unbetted_meanurf = [procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep tmp_unbetted_meanurf(1).name];
    [~,tmp_msg] = system([fsl_bin_path 'bet ' tmp_unbetted_meanurf ...
        ' ' regexprep(tmp_unbetted_meanurf,'mean','mmean') ...
        ' -f ' fsl_bet_default_th_avgvol ' -R -m']); % use -ROBUST
    if ~isempty(tmp_msg), disp(tmp_msg); end
    clear tmp_msg tmp_unbetted_meanurf;
    fprintf(['\nDone! (' num2str(toc) ' s)\n']);
else
    fprintf(['\nSkipping mean*.nii volume betting (already exists for ' which_visit_as_r_u_template ') for ' subj_ID '...\n']);
end
clear tmp_run_name;




%   5) Run PhysIO TAPAS in non-verbose mode in order to create the
%       regressors of interest related to motion (12, including
%       derivatives), physiological noise (16, from cardiac & breathing
%       belt), and FDs as specified from input (FD_th).
%   - Note that FDs are not regressed out but only kept for later usage;

all_NF_visits = dir([procpath subj_ID filesep NF_dir filesep 'V*']);
for v = 1:length(all_data_dirs)
    v_str = all_data_dirs(v).name((regexp(all_data_dirs(v).name,'GNV_V')+4):(regexp(all_data_dirs(v).name,'GNV_V')+6)); % find v_str from file name instead of assuming from v value
    all_NF_runs_per_v_str = dir([procpath subj_ID filesep NF_dir filesep v_str filesep 'Run*']);
    for r = 1:length(all_NF_runs_per_v_str)
        if ~exist([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep physiological_folder_name filesep],'dir')
            mkdir([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep physiological_folder_name filesep]);
        elseif exist([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep physiological_folder_name filesep 'multiple_regressors.txt'],'file')
            fprintf(['\nSkipping PhysIO TAPAS toolbox for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name ' (multiple_regressors.txt already exists).\n']);
            continue;
        end
        % manually addressing physiological recordings issues...
        if strcmp(subj_ID,'003GNV')
            if strcmp(v_str,'V11') && strcmp(all_NF_runs_per_v_str(r).name,'Run4_Transfer'), continue, end % corrupted fMRI data -> no need to run PhysIO
        elseif strcmp(subj_ID,'005GNV')
            if strcmp(v_str,'V06'), continue, end % whole visit with unstable pulse recordings unfortunately
        elseif strcmp(subj_ID,'008GNV')
            if strcmp(v_str,'V02') && (strcmp(all_NF_runs_per_v_str(r).name,'Run3') || strcmp(all_NF_runs_per_v_str(r).name,'Run5') || strcmp(all_NF_runs_per_v_str(r).name,'Run6') || strcmp(all_NF_runs_per_v_str(r).name,'Run7')), continue, end % problematic artefacts too (unstable pulse)
            if strcmp(v_str,'V06'), continue, end % whole visit with unstable pulse recordings unfortunately
            if strcmp(v_str,'V07') && strcmp(all_NF_runs_per_v_str(r).name,'Run2'), continue, end % unstable pulse
        elseif strcmp(subj_ID,'012GNV')
            if strcmp(v_str,'V02') && strcmp(all_NF_runs_per_v_str(r).name,'Run7'), continue, end % physio too artefacted
            if strcmp(v_str,'V04') && strcmp(all_NF_runs_per_v_str(r).name,'Run7'), continue, end % physio too artefacted
        elseif strcmp(subj_ID,'016GNV')
            %if strcmp(v_str,'V09') && strcmp(all_NF_runs_per_v_str(r).name,'Run2'), continue, end % cut
            %if strcmp(v_str,'V10') && strcmp(all_NF_runs_per_v_str(r).name,'Run5_Transfer'), continue, end % cut
            if strcmp(v_str,'V12') && (strcmp(all_NF_runs_per_v_str(r).name,'Run6') || strcmp(all_NF_runs_per_v_str(r).name,'Run7')), continue, end % cut
        elseif strcmp(subj_ID,'020GNV')
            if strcmp(v_str,'V03') && strcmp(all_NF_runs_per_v_str(r).name,'Run1'), continue, end % cut as well
        elseif strcmp(subj_ID,'021GNV')
            if strcmp(v_str,'V06') && strcmp(all_NF_runs_per_v_str(r).name,'Run6'), continue, end % physio too artefacted
            if strcmp(v_str,'V10') && strcmp(all_NF_runs_per_v_str(r).name,'Run1'), continue, end % not enough data (258/270 triggers)
            if strcmp(v_str,'V12') && strcmp(all_NF_runs_per_v_str(r).name,'Run7'), continue, end % physio too artefacted
            if strcmp(v_str,'V14') && (strcmp(all_NF_runs_per_v_str(r).name,'Run1') || strcmp(all_NF_runs_per_v_str(r).name,'Run4_Transfer') || strcmp(all_NF_runs_per_v_str(r).name,'Run7_NoFeedback')), continue, end % physio too artefacted or no data
            if strcmp(v_str,'V15') && (strcmp(all_NF_runs_per_v_str(r).name,'Run2') || strcmp(all_NF_runs_per_v_str(r).name,'Run3_Transfer') ...
                    || strcmp(all_NF_runs_per_v_str(r).name,'Run5_NoFeedback') || strcmp(all_NF_runs_per_v_str(r).name,'Run6_NoFeedback')), continue, end % physio too artefacted or no data
        elseif strcmp(subj_ID,'023GNV')
            if strcmp(v_str,'V13') && strcmp(all_NF_runs_per_v_str(r).name,'Run1'), continue, end % physio too artefacted
        elseif strcmp(subj_ID,'025GNV')
            if strcmp(v_str,'V15') && strcmp(all_NF_runs_per_v_str(r).name,'Run7_NoFeedback'), continue, end % physio too artefacted
        elseif strcmp(subj_ID,'026GNV')
            if strcmp(v_str,'V11') && strcmp(all_NF_runs_per_v_str(r).name,'Run6'), continue, end % physio too artefacted
        elseif strcmp(subj_ID,'034GNV')
            if strcmp(v_str,'V12') && strcmp(all_NF_runs_per_v_str(r).name,'Run5_Transfer'), continue, end % physio too artefacted
        elseif strcmp(subj_ID,'041GNV')
            if strcmp(v_str,'V04') && strcmp(all_NF_runs_per_v_str(r).name,'Run6'), continue, end % physio artefacted near the end...
        end
        %
        tmp_resp_pulse_mat_NF_dir = dir([basepath_data all_data_dirs(v).name filesep '*Aud_NF_Run_' all_NF_runs_per_v_str(r).name(4) '*']);
        if isempty(tmp_resp_pulse_mat_NF_dir)
            error(['Unable to find NF raw data directory corresponding to ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name '. Please check manually...']);
        end
        %tmp_resp_pulse_mat_file2 = dir([basepath_data all_data_dirs(v).name filesep ...
        %    tmp_resp_pulse_mat_NF_dir(1).name filesep all_data_dirs(v).name(1:17) '*_resp_pulse_Aud_NF_run' all_NF_runs_per_v_str(r).name(4) '.mat']); % here added (1:17) because sometimes dir names extend more than XXYYZZ_000GNV_V00
        tmp_resp_pulse_mat_file2 = dir([basepath_data all_data_dirs(v).name filesep ...
            tmp_resp_pulse_mat_NF_dir(1).name filesep all_data_dirs(v).name(1:17) '_resp_pulse_Aud_NF_run' all_NF_runs_per_v_str(r).name(4:end) '.mat']); % here added (1:17) because sometimes dir names extend more than XXYYZZ_000GNV_V00
        if isempty(tmp_resp_pulse_mat_file2)
            warning(['Unable to find .mat physiological recordings file for Run' all_NF_runs_per_v_str(r).name(4:end) ' in ' all_data_dirs(v).name filesep ...
                tmp_resp_pulse_mat_NF_dir(1).name '. Please check manually...']); % changed to warning; and then just skip the PhysIO part for this visit...
            continue; % skip this run
        end
        tmp_resp_pulse_mat_file = [basepath_data all_data_dirs(v).name filesep tmp_resp_pulse_mat_NF_dir(1).name filesep tmp_resp_pulse_mat_file2(1).name];
        tmp_physio_struct = load(tmp_resp_pulse_mat_file);
        if ~isequal(length(find(diff(tmp_physio_struct.data(:,3))))/2,270)
            fprintf(['\nWarning: Physio data for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name ' does not contain 270 scan peaks (but ' ...
                num2str(length(find(tmp_physio_struct.data(:,3)))/2) ')! Check manually...\n']);
        end
        if isequal(length(find(tmp_physio_struct.data(:,3)))/2,0) % take into account here the rare case in which triggers are missing -> either fix the .mat or skip for the moment!
            fprintf(['\nWarning: Triggers seem to be missing for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name ' (found ' ...
                num2str(length(find(tmp_physio_struct.data(:,3)))/2) ' peaks)! Either fix the .mat file manually and re-run script; or ignore physiological signals' ...
                ' regression for the moment (this option is chosen for now)...\n']);
            clear tmp_physio_struct tmp_resp_pulse_mat_file tmp_resp_pulse_mat_file2 tmp_resp_pulse_mat_NF_dir;
            continue; % skip to next run
        end
        fprintf(['\nFound ' num2str(length(find(tmp_physio_struct.data(:,3)))/2) ' peaks for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name '.\n']);
        if ((find(tmp_physio_struct.data(:,3),1,'last')-find(tmp_physio_struct.data(:,3),1,'first'))/6000)<6.725
            fprintf(['\nWarning: Physio data for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name ' seems to be less than 6.725 min! Check manually...\n']);
        end
        tmp_rp_movement_file = dir([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep 'rp_' aflag_prefix 'f*.txt']);
        if isempty(tmp_rp_movement_file)
            error(['Unable to find rp_*.txt file for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name '! Can''t proceed with motion regressors and FDs...']);
        end
        tmp_rp_movement_file = [procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep tmp_rp_movement_file(1).name];
        matlabbatch{1}.spm.tools.physio.save_dir = {[procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep physiological_folder_name]};
        matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Biopac_Mat';
        matlabbatch{1}.spm.tools.physio.log_files.cardiac = {tmp_resp_pulse_mat_file};
        matlabbatch{1}.spm.tools.physio.log_files.respiration = {tmp_resp_pulse_mat_file};
        matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
        matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
        matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = find(tmp_physio_struct.data(:,3),1)/100; % this needs to be computed as such (delay for starting the Biopac...)
        matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'first';
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = 64; % same as below (64 z-slices)
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = 1.5; % 1500 ms
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = 270; %length(dir([procpath subj_ID filesep NF_dir filesep v_str filesep all_NF_runs_per_v_str(r).name filesep 'f*.nii'])); % 270
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = 32; % selected half-z-slicing as onset slice
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
        matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
        %matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
        matlabbatch{1}.spm.tools.physio.scan_timing.sync.scan_timing_log = struct([]);
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
        matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
        matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt'; % output file
        matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat'; % same
        matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none'; % default
        matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = true;
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
        matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
        matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
        matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
        matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
        matlabbatch{1}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {tmp_rp_movement_file};
        matlabbatch{1}.spm.tools.physio.model.movement.yes.order = 12; % 6 motion regressors and their derivatives
        matlabbatch{1}.spm.tools.physio.model.movement.yes.censoring_method = 'FD';
        matlabbatch{1}.spm.tools.physio.model.movement.yes.censoring_threshold = FD_th; % compute FDs with this threshold
        matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
        matlabbatch{1}.spm.tools.physio.verbose.level = 0;
        matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = 'physioreg.fig'; % base name for figures
        matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
        fprintf(['\nRunning the TAPAS PhysIO preprocessing batch for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name '...\n']); tic;
        %spm('defaults','fmri');
        %spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        fprintf(['\n-> TAPAS PhysIO preprocessing for ' subj_ID ' ' v_str ' ' all_NF_runs_per_v_str(r).name ' done ! (' num2str(toc) ' s)\n']);
        clear matlabbatch tmp_physio_struct tmp_resp_pulse_mat_file* tmp_rp_movement_file;
    end
end




%   6) T1_V01 segmentation should be available from rsfMRI! Only running
%       CSF/WM masks remapping for further regression in the following step
%       (i.e. maps back the CSF & WM ROIs to functional space according to
%       the template obtained at step 1));

T1_V01_path = [procpath subj_ID filesep 'T1_V01' filesep];
if ~exist(T1_V01_path,'dir') || ~exist([T1_V01_path 'Segmented' filesep],'dir')
    error(['Unable to continue CSF and WM masks remapping for ' subj_ID ' (NF) because the T1_V01 folder or its segmentation does not exist!']);
end
tmp_T1_fname = dir([T1_V01_path 's*.nii']);
tmp_T1_fname = tmp_T1_fname(1).name; % here there shouldn't be any problems...
if urflag==1, if aflag, tmp_suffix = 'ua'; else, tmp_suffix = 'u'; end, elseif urflag==0, if aflag, tmp_suffix = 'a'; else, tmp_suffix = ''; end, end
tmp_run_name = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*']);
if isempty(tmp_run_name)
    error(['Unable to find the preprocessed Run1 folder for ' subj_ID ' ' which_visit_as_r_u_template '; maybe it does not exist for this visit? Please check manually...']);
end
tmp_run_name = tmp_run_name(1).name;
% convert the anatomical AAL90 to functional mean (NF) template
if isempty(dir([T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'c1*_AtlasAAL90_correctLR_th04_functional_mean' tmp_suffix 'f_NF.nii']))
    tmp_meanurfVXX = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean' tmp_suffix 'f*.nii']);
    if isempty(tmp_meanurfVXX), error(['Unable to find mean functional (pre-processed) volume for ' subj_ID ' ' ...
            which_visit_as_r_u_template ' ' tmp_run_name ', with pre-proc prefix ' r_u_prefix aflag_prefix '. Please check manually...']); end
    tmp_meanurfVXX = [procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep tmp_meanurfVXX(1).name];
    convert_atlas(...
        tmp_meanurfVXX,...
        [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'c1' tmp_T1_fname(1:end-4) '_AtlasAAL90_correctLR_th04.nii'],...
        [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'c1' tmp_T1_fname(1:end-4) '_AtlasAAL90_correctLR_th04_functional_mean' tmp_suffix 'f_NF.nii'],...
        'B'); % [not 'A' because of convert_atlas internals...]
else, fprintf(['\nAtlased functional c1 (to AAL90_correctLR, NF) already exists for ' subj_ID ' with th = 0.4! Skipping...\n']);
end
% reslicing c2 and c3 wrt the functional (atlased) c1 above
if ~exist([T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'rc2' tmp_T1_fname(1:end-4) '_NF.nii'],'file') || ...
        ~exist([T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'rc3' tmp_T1_fname(1:end-4) '_NF.nii'],'file')
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.realign.write.data{1,1} = ...
        [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'c1' tmp_T1_fname(1:end-4) '_AtlasAAL90_correctLR_th04_functional_mean' tmp_suffix 'f_NF.nii,1'];
    matlabbatch{1}.spm.spatial.realign.write.data{2,1} = ...
        [T1_V01_path 'Segmented' filesep 'c2' tmp_T1_fname ',1'];
    matlabbatch{1}.spm.spatial.realign.write.data{3,1} = ...
        [T1_V01_path 'Segmented' filesep 'c3' tmp_T1_fname ',1'];
    matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
    %spm('defaults','fmri'); spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    movefile([T1_V01_path 'Segmented' filesep matlabbatch{1}.spm.spatial.realign.write.roptions.prefix 'c2' tmp_T1_fname],...
        [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep matlabbatch{1}.spm.spatial.realign.write.roptions.prefix 'c2' tmp_T1_fname(1:end-4) '_NF.nii']);
    movefile([T1_V01_path 'Segmented' filesep matlabbatch{1}.spm.spatial.realign.write.roptions.prefix 'c3' tmp_T1_fname],...
        [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep matlabbatch{1}.spm.spatial.realign.write.roptions.prefix 'c3' tmp_T1_fname(1:end-4) '_NF.nii']);
else, fprintf(['\nResliced c2 and c3 (NF) already exist for ' subj_ID '! Skipping...\n']);
end
% GM & WM masks deformation & mapping
if ~exist([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep],'dir'), mkdir([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep]); end
if ~exist([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'rw' WM_template_name_func],'file') || ...
        ~exist([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'rw' CSF_template_name_func],'file') || ...
        ~exist([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'irwc2' WM_template_name_func],'file') || ...
        ~exist([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'irwc3' CSF_template_name_func],'file')
        %~exist([T1_V01_path 'Segmented' filesep 'irwc2' WM_template_name_func],'file') || ...
        %~exist([T1_V01_path 'Segmented' filesep 'irwc3' CSF_template_name_func],'file')
    clear matlabbatch;
    if ~exist([T1_V01_path 'Segmented' filesep 'y_' tmp_T1_fname],'file')
        error(['Unable to find the y_*.nii forward transform (T1 V01) for participant ' subj_ID '... Please recheck the entire routine!']);
    end
    tmp_meanurfVXX = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean' tmp_suffix 'f*.nii']);
    if isempty(tmp_meanurfVXX), error(['Unable to find mean functional (pre-processed) volume for ' subj_ID ' ' ...
            which_visit_as_r_u_template ' ' tmp_run_name ', with pre-proc prefix ' r_u_prefix aflag_prefix '. Please check manually...']); end
    tmp_meanurfVXX = [procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep tmp_meanurfVXX(1).name];
    matlabbatch{1}.spm.util.defs.comp{1}.def = {[T1_V01_path 'Segmented' filesep 'y_' tmp_T1_fname]};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1,1} = [templates_wmcsf_path CSF_template_name];
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames{2,1} = [templates_wmcsf_path WM_template_name];
    matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {[T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {tmp_meanurfVXX};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix = 'rw';
    %spm('defaults','fmri'); spm_jobman('initcfg');
    spm_jobman('run',matlabbatch); % it will deform the templates CSF & WM masks and reslice them to functional resolution
    tmp_outprefix = matlabbatch{1}.spm.util.defs.out{1}.push.prefix;
    movefile([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix WM_template_name],...
        [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix WM_template_name_func]);
    movefile([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix CSF_template_name],...
        [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix CSF_template_name_func]);
    clear matlabbatch;
    % WM masking with rc2* segmentation from before (output : irwc2White*)
    matlabbatch{1}.spm.util.imcalc.input{1,1} = [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix WM_template_name_func ',1'];
    matlabbatch{1}.spm.util.imcalc.input{2,1} = [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'rc2' tmp_T1_fname(1:end-4) '_NF.nii,1'];
    matlabbatch{1}.spm.util.imcalc.output = '';
    matlabbatch{1}.spm.util.imcalc.outdir = {[T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF']};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*(i2>0)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -4;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    % CSF masking with rc3* segmentation from before (output : irwc3Csf*)
    matlabbatch{2}.spm.util.imcalc.input{1,1} = [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep tmp_outprefix CSF_template_name_func ',1'];
    matlabbatch{2}.spm.util.imcalc.input{2,1} = [T1_V01_path 'Atlased_AAL90_NewAutoLabellingOut' filesep 'rc3' tmp_T1_fname(1:end-4) '_NF.nii,1'];
    matlabbatch{2}.spm.util.imcalc.output = '';
    matlabbatch{2}.spm.util.imcalc.outdir = {[T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF']};
    matlabbatch{2}.spm.util.imcalc.expression = 'i1.*(i2>0)';
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = -4;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    %spm('defaults','fmri'); spm_jobman('initcfg');
    spm_jobman('run',matlabbatch); % it will mask the previous functional masks with resp. rc2 and rc3 (functional) segmentations
    movefile([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'i' tmp_outprefix WM_template_name_func],...
        [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'i' tmp_outprefix 'c2' WM_template_name_func]);
    movefile([T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'i' tmp_outprefix CSF_template_name_func],...
        [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep 'i' tmp_outprefix 'c3' CSF_template_name_func]);
    %clear tmp_outprefix matlabbatch tmp_meanurfV01 tmp_suffix;
else, fprintf(['\nMasked and resliced (NF) irwc2White* and irwc3Csf* masks already exist for ' subj_ID '! Skipping...\n']);
end
clear tmp_outprefix	matlabbatch tmp_meanurfV01 tmp_T1_fname; % don't clear tmp_run_name yet! used later
%
% % check if parc_fs exists from DWI/FreeSurfer extraction, so that we also
% %     % reslice the high-res parc_fs to parc_fs_functional
% %     if exist([procpath subj_ID filesep 'DTI' filesep 'V01' filesep 'parc_fs.mif'],'file')
% %         if ~exist([T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep],'dir')
% %             mkdir([T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep]);
% %         end
% %         if ~exist([T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep 'parc_fs_' subj_ID '_V01.nii'],'file')
% %             [~,tmpmsg] = system([filesep 'usr' filesep 'local' filesep 'bin' filesep 'mrconvert ' procpath subj_ID filesep 'DTI' ...
% %                 filesep 'V01' filesep 'parc_fs.mif' ' ' T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep 'parc_fs_' subj_ID '_V01.nii']);
% %             if ~isempty(tmpmsg), disp(tmpmsg); end
% %         end
% %         if ~exist([T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep 'parc_fs_' subj_ID '_V01_' num2str(func_vols_size(1)) 'x' ...
% %                 num2str(func_vols_size(2)) 'x' num2str(func_vols_size(3)) '.nii'],'file')
% %             tmp_meanunewprefVXX = dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep 'meanunewpref*.nii']);
% %             if isempty(tmp_meanunewprefVXX), error(['Unable to find mean functional (pre-proc) volume for ' subj_ID ' V01 Run1 (with pre-proc prefix ' proc_files_prefix '?)!']); end
% %             tmp_meanunewprefVXX = [procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep tmp_meanunewprefVXX(1).name];
% %             convert_atlas(tmp_meanunewprefVXX,...
% %                 [T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep 'parc_fs_' subj_ID '_V01.nii'],...
% %                 [T1_VXX_path 'fs_parc_to_functional_mean' proc_files_prefix(1:end-1) '_NF' filesep 'parc_fs_' subj_ID '_V01_' num2str(func_vols_size(1)) 'x' ...
% %                 num2str(func_vols_size(2)) 'x' num2str(func_vols_size(3)) '.nii'],...
% %                 'B'); % temp, just because not AAL
% %         else, fprintf(['\nparc_fs_' subj_ID '_V01_' num2str(func_vols_size(1)) 'x' ...
% %                 num2str(func_vols_size(2)) 'x' num2str(func_vols_size(3)) '.nii already exists! Skipping...\n']);
% %         end
% %     else, fprintf(['\nUnable to find pars_fs.mif for ' subj_ID '... Skipping!\n']);
% %     end
%




%   7) Proper regression of nuisance variables as mentioned above;

if isempty(dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep reg_vols_prefix '*.nii'])) % no Reg2_*.nii (or reg_vols_prefix) files exist -> run the regression anew
    fprintf(['\nRunning regression on pre-processed volumes (NF) for ' subj_ID '...\n']); tmp_time = tic;
    CSF_mask_name = 'irwc3CsfMask_07_108x108x64.nii';
    WM_mask_name = 'irwc2WhiteMask_09_108x108x64.nii';
    irwCSF_path = [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep CSF_mask_name];
    irwWM_path = [T1_V01_path 'GM_CSF_to_functional_mean' tmp_suffix 'f_NF' filesep WM_mask_name];
    %
    fprintf(['\nLoading CSF mask... (' CSF_mask_name ')\n']);
    CSF_mask = load_untouch_nii(irwCSF_path);
    CSF_mask_data = double(CSF_mask.img); %single() is fine
    %
    fprintf(['\nLoading WM mask... (' WM_mask_name ')\n']);
    WM_mask = load_untouch_nii(irwWM_path);
    WM_mask_data = double(WM_mask.img); %single() is fine
    %
    try srv_name = char(java.net.InetAddress.getLocalHost.getHostName); catch, srv_name = 'unknown srv'; end
    try
        clear tmp_pcluster;
        tmp_pcluster = parcluster;
        tmp_pcluster.NumWorkers = num_workers_parallel;
        tmp_pool = parpool(num_workers_parallel);
        fprintf(['\nInitialized parallel pool on ' srv_name ' with ' num2str(tmp_pool.NumWorkers) ' workers...\n']);
    catch %, parpool; fprintf('\nDefault parpool initialized...\n'); % default
    end
    %
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            fprintf(['\nLoading NF time series... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']); tic;
            all_NF_run_volumes = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep r_u_prefix aflag_prefix 'f*.nii']);
            if ~isequal(length(all_NF_run_volumes),270)
                error(['Incorrect amount of NF pre-processed volumes found for ' subj_ID ' ' ...
                    all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ...
                    ' (found ' num2str(length(all_NF_run_volumes)) ' instead of 270). Please check manually!']);
            end
            all_NF_visit_data = zeros(func_vols_size(1),func_vols_size(2),func_vols_size(3),length(all_NF_run_volumes)); % 108x108x64x270
            for rr = 1:length(all_NF_run_volumes)
                all_NF_visit_data_vol = load_untouch_nii([procpath subj_ID filesep NF_dir filesep ...
                    all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_NF_run_volumes(rr).name]);
                all_NF_visit_data(:,:,:,rr) = double(all_NF_visit_data_vol.img);
            end, clear all_NF_visit_data_vol;
            %
            fprintf(['\nComputing average CSF and WM signals... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']);
            tmp_alldata = reshape(all_NF_visit_data,prod(func_vols_size),length(all_NF_run_volumes));
            CSF_avg = mean(tmp_alldata(find(CSF_mask_data(:)),:)); %#ok ;
            WM_avg = mean(tmp_alldata(find(WM_mask_data(:)),:)); %#ok ;
            if which_regression_approach==1
                CSF_avg_demean = CSF_avg - mean(CSF_avg);
                WM_avg_demean = WM_avg - mean(WM_avg);
            elseif which_regression_approach==2
                CSF_avg_diff = [0,diff(CSF_avg)];
                WM_avg_diff = [0,diff(WM_avg)]; % include derivatives of WM & CSF to be regressed out in 2nd approach
            end, clear tmp_alldata;
            %   building the regression matrix XREG (timepoints x regressors); so in
            %   principle, e.g. 370x37 if having 3 trends (cst, lin, qdrtic),
            %   CSF & WM avg and demeaned and their derivatives, 18 physiological
            %   signals obtained using RETROICOR (PhysIO TAPAS Toolbox for SPM),
            %   12 motion parameters (rp) and their derivatives...
            fprintf(['\nBuilding the regressors... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']);
            if ~exist([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                    all_NF_runs_per_visit(r).name filesep physiological_folder_name filesep 'multiple_regressors.txt'],'file')
                fprintf(['\nWARNING: Multiple regressors file NOT FOUND for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '!\nSkipping PhysIO regressors (only 12 motion regressors will be used)!\n']);
                tmp_additional_reg_data = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep 'rp_' aflag_prefix 'f*.txt']);
                if isempty(tmp_additional_reg_data), error(['\nrp*.txt also NOT FOUND for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '! Aborting...']); end
                tmp_additional_reg_data = load([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep tmp_additional_reg_data(1).name]);
                outvols_data_prefix = [reg_vols_prefix 'MotionOnly_']; % reg_vols_prefix is usually 'Reg2_' or so
            else % physiological regressors available
                tmp_additional_reg_data = load([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep physiological_folder_name filesep 'multiple_regressors.txt']);
                outvols_data_prefix = [reg_vols_prefix 'Physio_'];
            end
            clear XREG;
            if which_regression_approach==1
                XREG = [ones(length(all_NF_run_volumes),1), (1:length(all_NF_run_volumes))'/length(all_NF_run_volumes), ...
                    ((1:length(all_NF_run_volumes))'/length(all_NF_run_volumes)).^2, WM_avg_demean', CSF_avg_demean', tmp_additional_reg_data];
            elseif which_regression_approach==2 % in this approach we have 30 regressors in total to add (18 physiological : 6+8+4 interaction; 12 motion (6 + 6 derivatives))
                % but the above is only true if the cardiac regressors could be created
                if size(tmp_additional_reg_data,2)>=30 % assuming here that cardiac regressors could be created
                    tmp_additional_reg_data = tmp_additional_reg_data(:,1:30); % /!\ shrink to remove FDs regressors /!\
                elseif regexp(outvols_data_prefix,'MotionOnly','once') % in this case, the file multiple_regressors.txt does not exist and therefore we only have rp*.txt as 270x6 motion parameters
                    tmp_additional_reg_data = [tmp_additional_reg_data diff([zeros(1,6);tmp_additional_reg_data])]; %#ok % make it 270x12 with the derivatives
                else
                    error(['Unable to properly evaluate ''tmp_additional_reg_data'' variable''s length for ' ...
                        subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '. Double-check manually!']); % should not happen
                end
                XREG = [ones(length(all_NF_run_volumes),1), (1:length(all_NF_run_volumes))'/length(all_NF_run_volumes), ...
                    ((1:length(all_NF_run_volumes))'/length(all_NF_run_volumes)).^2, WM_avg', CSF_avg', WM_avg_diff', CSF_avg_diff', tmp_additional_reg_data];
            end
            %
            fprintf(['\nRegressing out... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ', reg. approach #' num2str(which_regression_approach) ')\n']);
            all_NF_visit_data_reg = zeros(func_vols_size(1),func_vols_size(2),func_vols_size(3),length(all_NF_run_volumes)); % 108x108x64x270
            parfor i = 1:func_vols_size(1)
                for j = 1:108%func_vols_size(2) % [parfor]
                    for k = 1:64%func_vols_size(3) % [parfor]
                        %[beta,res,SSE,SSR,T] =
                        %y_regress_ss(squeeze(V0(i,j,k,V0idx)),X);
                        [~,residuals] = y_regress_ss(squeeze(all_NF_visit_data(i,j,k,:)),XREG);
                        all_NF_visit_data_reg(i,j,k,:) = residuals + mean(squeeze(all_NF_visit_data(i,j,k,:))); % put back the mean signal
                    end
                end
            end
            if any(isnan(all_NF_visit_data_reg)) % should not happen as well but for sanity
                fprintf(['\nFound NaNs!! Removing NaNs from regressed-out matrix... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']);
                all_NF_visit_data_reg(isnan(all_NF_visit_data_reg)) = 0;
            end
            %
            fprintf(['\nMasking and writing regressed-out data... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ', with prefix ''' outvols_data_prefix ''')\n']);
            tmp_masking_vol = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mmean*_mask.nii']);
            if isempty(tmp_masking_vol), error(['Unable to find average EPI volume mask (mmean*_mask.nii) for ' subj_ID ' ' which_visit_as_r_u_template ' ' tmp_run_name '... Please check manually!']); end
            tmp_masking_vol = load_untouch_nii([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep tmp_masking_vol(1).name]);
            tmp_masking_vol_data = logical(tmp_masking_vol.img);
            all_NF_visit_data_reg = reshape(all_NF_visit_data_reg,prod(func_vols_size),length(all_NF_run_volumes));
            all_NF_visit_data_reg(find(tmp_masking_vol_data==0),:) = 0; %#ok
            all_NF_visit_data_reg = reshape(all_NF_visit_data_reg,func_vols_size(1),func_vols_size(2),func_vols_size(3),length(all_NF_run_volumes));
            for rr = 1:length(all_NF_run_volumes)
                all_NF_visit_data_vol = load_untouch_nii([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_NF_run_volumes(rr).name]);
                all_NF_visit_data_vol.img = uint16(all_NF_visit_data_reg(:,:,:,rr)); % cast back the data to uint16
                all_NF_visit_data_vol.fileprefix = [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep ...
                    outvols_data_prefix all_NF_run_volumes(rr).name(1:end-4)]; % [(1:end-4) because we remove the .nii extension for the fileprefix]
                save_untouch_nii(all_NF_visit_data_vol,[all_NF_visit_data_vol.fileprefix '.nii']);
            end, clear all_NF_visit_data_vol;
            fprintf(['\nDone for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '... (in ' num2str(toc) ' s)\n']);
        end
    end
    fprintf(['\nRegression (NF) done for ' subj_ID '... (in ' num2str(toc(tmp_time)) ' s)\n']); clear tmp_time;
else
    fprintf(['\nRegressed volumes (NF, with prefix ' reg_vols_prefix ') seem to already exist for ' subj_ID '! Skipping...\n']);
end




%   8) ROIs unwarping/realignment (according to a specific visit);
%       This is useful for later on processing/analysis;

remapped_ROIs_NF_dir = [procpath subj_ID filesep NF_dir filesep 'remapped_ROIs' filesep];
if ~exist(remapped_ROIs_NF_dir,'dir'), mkdir(remapped_ROIs_NF_dir); end
if isempty(dir([remapped_ROIs_NF_dir 'mr*ROI_*.nii'])) % then no realigned/unwarped ROIs are to be found, proceed with this part
    fprintf(['\nWorking on NF ROIs (unwarping/realignment) for ' subj_ID '...\n']);
    tmp_defaultROIsvisit = 'V01';
    % special cases here >
    if strcmp(subj_ID,'003GNV'), tmp_defaultROIsvisit = 'V02'; end
    if strcmp(subj_ID,'005GNV'), tmp_defaultROIsvisit = 'V03'; end
    if strcmp(subj_ID,'016GNV'), tmp_defaultROIsvisit = 'V02'; end
    if strcmp(subj_ID,'021GNV'), tmp_defaultROIsvisit = 'V03'; end % because used 002GNV's ROIs for V01, V02 and V05 (attempt to come back to larger ROIs)
    if strcmp(subj_ID,'034GNV'), tmp_defaultROIsvisit = 'V02'; end
    % <
    tmp_AudROIs = dir([basepath_data '*' subj_ID '_' tmp_defaultROIsvisit filesep '*Auditorydownreg_ROIs_1*' filesep '*ROI*.nii']);
    if ~isequal(length(tmp_AudROIs),3), error(['Unable to find all the NF ROIs (' tmp_defaultROIsvisit ') for ' subj_ID '... Please check manually!']); end
    tmp_EPItemplate = dir([basepath_data '*' subj_ID '_' tmp_defaultROIsvisit filesep '*Auditorydownreg_EPI_Template_1*' filesep 'meanf*.nii']);
    if isempty(tmp_EPItemplate), error(['Unable to find NF (experimental) EPI template (' tmp_defaultROIsvisit ') for ' subj_ID '! Please check manually!']); end
    %tmp_EPItemplate = [tmp_EPItemplate(1).folder filesep tmp_EPItemplate(1).name];
    for r = 1:length(tmp_AudROIs), copyfile([tmp_AudROIs(r).folder filesep tmp_AudROIs(r).name],[remapped_ROIs_NF_dir tmp_AudROIs(r).name]); end
    copyfile([tmp_EPItemplate(1).folder filesep tmp_EPItemplate(1).name],[remapped_ROIs_NF_dir tmp_EPItemplate(1).name]);
    if urflag==1
        tmp_EPIref = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*' filesep 'mean' r_u_prefix aflag_prefix '*.nii']);
        if isempty(tmp_EPIref), error(['Unable to find NF preprocessed (mean' r_u_prefix aflag_prefix '*.nii) template for ' subj_ID ' ' which_visit_as_r_u_template '. Please double-check!']); end
        tmp_EPIref = [tmp_EPIref(1).folder filesep tmp_EPIref(1).name];
        tmp_vdm_file = dir([basepath_data '*' subj_ID '_' tmp_defaultROIsvisit filesep '_posthoc' filesep 'gre_field_nii' filesep 'vdm5_*.nii']);
        if isempty(tmp_vdm_file), error(['Unable to find vdm5 (gre_field .nii) file (' tmp_defaultROIsvisit ') for ' subj_ID '... Please check manually!']); end
        tmp_vdm_file = [tmp_vdm_file(1).folder filesep tmp_vdm_file(1).name];
        clear matlabbatch;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.scans = ...
            {
            [remapped_ROIs_NF_dir tmp_EPItemplate(1).name ',1']
            [remapped_ROIs_NF_dir tmp_AudROIs(1).name ',1']
            [remapped_ROIs_NF_dir tmp_AudROIs(2).name ',1']
            [remapped_ROIs_NF_dir tmp_AudROIs(3).name ',1']
            };
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.vdmfile{1} = [tmp_vdm_file ',1'];
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2; % posterior-anterior dist. direction
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [2 0]; % don't create avg image from the ROIs
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.mask = 0; % for now, don't mask
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix = r_u_prefix; % 'u'
        matlabbatch{2}.spm.spatial.coreg.estwrite.ref = {[tmp_EPIref ',1']};
        matlabbatch{2}.spm.spatial.coreg.estwrite.source = {[remapped_ROIs_NF_dir ...
            matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix tmp_EPItemplate(1).name ',1']};
        matlabbatch{2}.spm.spatial.coreg.estwrite.other = {
            [remapped_ROIs_NF_dir matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix tmp_AudROIs(1).name ',1']
            [remapped_ROIs_NF_dir matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix tmp_AudROIs(2).name ',1']
            [remapped_ROIs_NF_dir matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix tmp_AudROIs(3).name ',1']
            };
        matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ncc'; % try norm. cross corr.
        matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        fprintf(['Unwarping ' subj_ID ' ' tmp_defaultROIsvisit ' (experimental) mean EPI template & ROIs...\n']); tic;
        spm_jobman('run',matlabbatch);
        fprintf(['Done (' num2str(toc) ' s).\nMasking the unwarped ROIs (' subj_ID ' ' tmp_defaultROIsvisit ') with a ' reg_vols_prefix '* volume...\n']);
        tmp_maskingprefix = [matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix];
    else
        tmp_EPIref = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*' filesep 'mean' aflag_prefix '*.nii']);
        if isempty(tmp_EPIref), error(['Unable to find NF preprocessed (mean' aflag_prefix '*.nii) template for ' subj_ID ' ' which_visit_as_r_u_template '. Please double-check!']); end
        tmp_EPIref = [tmp_EPIref(1).folder filesep tmp_EPIref(1).name];
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[tmp_EPIref ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[remapped_ROIs_NF_dir tmp_EPItemplate(1).name ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
            [remapped_ROIs_NF_dir tmp_AudROIs(1).name ',1']
            [remapped_ROIs_NF_dir tmp_AudROIs(2).name ',1']
            [remapped_ROIs_NF_dir tmp_AudROIs(3).name ',1']
            };
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ncc';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = r_u_prefix; % 'r'
        fprintf(['Realigning ' subj_ID ' ' tmp_defaultROIsvisit ' (experimental) mean EPI template & ROIs...\n']); tic;
        spm_jobman('run',matlabbatch);
        fprintf(['Done (' num2str(toc) ' s).\nMasking the realigned ROIs (' subj_ID ' ' tmp_defaultROIsvisit ') with a ' reg_vols_prefix '* volume...\n']);
        tmp_maskingprefix = r_u_prefix;
    end
    tmp_first_template = dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep reg_vols_prefix '*.nii']); % here it's fine to let V01 and Run1
    if isempty(tmp_first_template), error(['Unable to find any ' reg_vols_prefix '* files (NF) for ' subj_ID '! Please check manually...']); end
    tmp_first_template = load_untouch_nii([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep tmp_first_template(1).name]); % same
    tmp_first_template = single(tmp_first_template.img(:));
    for r = 1:length(tmp_AudROIs)
        tmp_ROI = load_untouch_nii([remapped_ROIs_NF_dir tmp_maskingprefix tmp_AudROIs(r).name]); % prefixed unwarped ROIs to load
        tmp_ROI_data = single(tmp_ROI.img(:));
        if length(unique(tmp_ROI_data))~=2
            fprintf(['Correcting ' tmp_maskingprefix tmp_AudROIs(r).name ' for un-binary values...\n']);
            tmp_ROI_data(tmp_ROI_data>1) = 1; tmp_ROI_data(tmp_ROI_data<0) = 0;
            tmp_ROI_data = round(tmp_ROI_data);
        end
        if ~isempty(find(isnan(tmp_ROI_data),1))
            fprintf(['Correcting ' tmp_maskingprefix tmp_AudROIs(r).name ' for NaN values (likely from unwarping)...\n']);
            tmp_ROI_data(isnan(tmp_ROI_data)) = 0;
        end
        tmp_ROI_data(tmp_first_template==0) = 0; % mask with first Reg* vol
        tmp_ROI_size = length(find(tmp_ROI_data>0));
        tmp_ROI.img = reshape(tmp_ROI_data,size(tmp_ROI.img,1),size(tmp_ROI.img,2),size(tmp_ROI.img,3));
        switch r
            case 1, tmp_ROI.fileprefix = [remapped_ROIs_NF_dir 'm' tmp_maskingprefix ...
                        'AudROI_L_' num2str(tmp_ROI_size)]; AudROI_L_ind = find(tmp_ROI_data>0);
            case 2, tmp_ROI.fileprefix = [remapped_ROIs_NF_dir 'm' tmp_maskingprefix ...
                        'AudROI_R_' num2str(tmp_ROI_size)]; AudROI_R_ind = find(tmp_ROI_data>0);
            case 3, tmp_ROI.fileprefix = [remapped_ROIs_NF_dir 'm' tmp_maskingprefix ...
                        'ctrlROI_' num2str(tmp_ROI_size)]; ctrlROI_ind = find(tmp_ROI_data>0);
        end
        save_untouch_nii(tmp_ROI,[tmp_ROI.fileprefix '.nii']);
        clear tmp_ROI tmp_ROI_size tmp_ROI_data;
        if urflag==1
            delete([remapped_ROIs_NF_dir matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix tmp_AudROIs(r).name]); % u*
        end
        delete([remapped_ROIs_NF_dir tmp_maskingprefix tmp_AudROIs(r).name]); % ru* or r*
    end
    save([remapped_ROIs_NF_dir 'ROIs_new_indices_' subj_ID '_NF.mat'],'AudROI_L_ind','AudROI_R_ind','ctrlROI_ind'); % save anyway
    clear matlabbatch tmp_first_template tmp_maskingprefix tmp_AudROIs tmp_EPIref tmp_EPItemplate tmp_vdm_file;
else
    fprintf(['\nPreprocessed ROIs seem to already exist (either unwarped or realigned) for ' subj_ID ' (NF)... Skipping...\n']);
end




%   9) Smoothing of Reg2_* volumes;

if sflag
    if isempty(dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(1).name ...
            filesep 'Run1' filesep 's' num2str(smoothing_factor) reg_vols_prefix '*.nii']))
        for v = 1:length(all_NF_visits)
            all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
            for r = 1:length(all_NF_runs_per_visit)
                all_NF_Reg_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                    all_NF_runs_per_visit(r).name filesep reg_vols_prefix '*.nii']);
                if ~isequal(length(all_NF_Reg_files),270)
                    error(['Incorrect number of ' reg_vols_prefix '*.nii files found in directory ' all_NF_visits(v).name ' ' ...
                        all_NF_runs_per_visit(r).name ' for ' subj_ID '! Found ' num2str(length(all_NF_Reg_files)) ' files instead of 270. Please check manually.']);
                end
                %if exist([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                %        all_NF_runs_per_visit(r).name filesep 's' num2str(smoothing_factor) all_NF_Reg_files(1).name],'file')
                %    fprintf(['\nSkipping smoothing for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ' (already exists)...\n']);
                %else
                clear matlabbatch;
                for rr = 1:length(all_NF_Reg_files)
                    matlabbatch{1}.spm.spatial.smooth.data{rr,1} = [procpath subj_ID filesep NF_dir filesep ...
                        all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_NF_Reg_files(rr).name ',1'];
                end
                matlabbatch{1}.spm.spatial.smooth.fwhm = [smoothing_factor smoothing_factor smoothing_factor];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0; % rather keep to 0
                matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(smoothing_factor)];
                %spm('defaults','fmri');
                %spm_jobman('initcfg');
                fprintf(['\nSmoothing... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']); tic;
                spm_jobman('run',matlabbatch);
                fprintf(['\nDone for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '... (in ' num2str(toc) ' s)\n']);
                %end
            end
        end
    else
        fprintf(['\nSmoothed volumes (NF, with prefix s' num2str(smoothing_factor) reg_vols_prefix ') seem to already exist for ' subj_ID '! Skipping...\n']);
    end
else
    fprintf(['\nSkipping smoothing (s' num2str(smoothing_factor) ' kernel) for ' subj_ID ' (sflag = 0).\n']); %#ok
end




%   9.1) Intersect mask creation for smoothed Reg2_* volumes (s6Reg2_*);

if ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_s' num2str(smoothing_factor) reg_vols_prefix(1:end-1) '.nii'],'file') % [end-1 because we want Reg2 and not Reg2_]
    fprintf(['\nRunning intersect mask creation (NF, for smoothed volumes (s' num2str(smoothing_factor) reg_vols_prefix(1:end-1) '*)) for ' subj_ID '...\n']); tic;
    all_NF_visits = dir([procpath subj_ID filesep NF_dir filesep 'V*']);
    tmp_run_name = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep 'Run1*']);
    if isempty(tmp_run_name)
        error(['Unable to find the preprocessed Run1 folder for ' subj_ID ' ' which_visit_as_r_u_template '; maybe it does not exist for this visit? Please check manually...']);
    end
    tmp_run_name = tmp_run_name(1).name;
    intersect_mask = dir([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep 'mean*.nii']); % does not matter whether meanuf/meanf
    if isempty(intersect_mask)
        error(['Unable to find (' r_u_prefix aflag_prefix ')mean EPI template for ' subj_ID ' (' which_visit_as_r_u_template ' ' ...
            tmp_run_name ') for intersect mask creation (for smoothed volumes)... Please double-check manually...']);
    end
    intersect_mask = load_untouch_nii([procpath subj_ID filesep NF_dir filesep which_visit_as_r_u_template filesep tmp_run_name filesep intersect_mask(1).name]);
    intersect_mask_data = double(intersect_mask.img); % [#ok]
    intersect_mask_data = intersect_mask_data(:);
    intersect_mask.fileprefix = ['intersect_mask_NF_' subj_ID '_s' num2str(smoothing_factor) reg_vols_prefix(1:end-1)];
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            all_s_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep 's' num2str(smoothing_factor) reg_vols_prefix '*.nii']);
            if ~isequal(length(all_s_files),270)
                error(['Incorrect number of s' num2str(smoothing_factor) reg_vols_prefix '*.nii volumes found in folder ' all_NF_visits(v).name ' ' ...
                    all_NF_runs_per_visit(r).name ' for ' subj_ID ' during (smoothed volumes) intersect mask creation. Please check manually!']);
            end
            for f = 1:length(all_s_files)
                tmp_vol_for_erosion = load_untouch_nii(...
                    [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_s_files(f).name]);
                tmp_vol_for_erosion = reshape(double(tmp_vol_for_erosion.img),size(tmp_vol_for_erosion.img,1)*...
                    size(tmp_vol_for_erosion.img,2)*size(tmp_vol_for_erosion.img,3),1);
                %intersect_mask_data = max(min(),0); % max precaution?
                intersect_mask_data = min(intersect_mask_data,tmp_vol_for_erosion);
            end
        end
    end
    intersect_mask.img = int16(reshape(intersect_mask_data,size(tmp_vol_for_erosion,1),size(tmp_vol_for_erosion,2),size(tmp_vol_for_erosion,3)));
    save_untouch_nii(intersect_mask,[procpath subj_ID filesep NF_dir filesep '_misc' filesep intersect_mask.fileprefix '.nii']);
    clear intersect_mask* tmp_vol_for_erosion all_s_files all_NF_runs_per_visit;
    fprintf(['\nDone (' subj_ID ', ' num2str(toc) ' s).\n']);
else
    fprintf(['\nSkipping smoothed volumes intersect mask creation (already exists) for ' subj_ID '...\n']);
end




%   10) Normalization to MNI space (of Reg2_* volumes);

if nflag
    fprintf(['\nRunning MNI normalization (NF) for ' subj_ID '...\n']); tic;
    tmp_T1_Vfolder = 'T1_V01';
    tmp_normalize_prefix = 'w';
    %if strcmp(subj_ID,'026GNV'), tmp_T1_Vfolder = 'T1_V02'; end
    %if strcmp(subj_ID,'037GNV'), tmp_T1_Vfolder = 'T1_V05'; end
    T1_pushforwardmap = dir([procpath subj_ID filesep tmp_T1_Vfolder filesep 'Segmented' filesep 'y_s*.nii']);
    if isempty(T1_pushforwardmap)
        error(['Unable to find T1 push-forward map (y_s*.nii) in the ''Segmented'' folder (' tmp_T1_Vfolder ') for ' subj_ID '! Aborting...']);
    else
        T1_pushforwardmap = [procpath subj_ID filesep tmp_T1_Vfolder filesep 'Segmented' filesep T1_pushforwardmap(1).name];
    end
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            if ~isempty(dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep ...
                    tmp_normalize_prefix reg_vols_prefix '*.nii'])) % e.g. wReg2_*.nii
                fprintf(['\nNormalized files (' tmp_normalize_prefix reg_vols_prefix '*.nii) seem to already exist for ' ...
                    subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name '... Skipping!\n']);
                continue;
            else
                fprintf(['\nNormalizing ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ' (' reg_vols_prefix '*.nii files) to MNI space...\n']);
                clear matlabbatch;
                all_NF_Reg_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name ...
                    filesep all_NF_runs_per_visit(r).name filesep reg_vols_prefix '*.nii']);
                if ~isequal(length(all_NF_Reg_files),270) % [do not allow 260 if skipping 5 at the beginning and 5 at the end]
                    error(['Incorrect number of ' reg_vols_prefix '*.nii files found in folder ' ...
                        all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep ' for ' subj_ID '! Please check manually.']);
                end
                matlabbatch{1}.spm.spatial.normalise.write.subj.def = {T1_pushforwardmap};
                matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cell(270,1);
                for f = 1:length(all_NF_Reg_files)
                    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{f,1} = [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                        all_NF_runs_per_visit(r).name filesep all_NF_Reg_files(f).name ',1'];
                end
                %matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
                matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72; 90 90 108];
                matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
                matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = tmp_normalize_prefix;
                %spm('defaults','fmri'); spm_jobman('initcfg');
                spm_jobman('run',matlabbatch);
            end
        end
    end
    clear tmp_T1_Vfolder matlabbatch;
    fprintf(['\nMNI normalization (NF, ' reg_vols_prefix '*.nii inputs) done for ' subj_ID ' (in ' num2str(toc) ' s).\n']);
else
    fprintf(['\nSkipping MNI normalization of ' reg_vols_prefix '*.nii volumes for ' subj_ID ' (nflag = 0).\n']); %#ok
end




%   10.1) Intersect mask creation for normalized Reg2_* volumes (wReg2_*);

if nflag
    if ~exist([procpath subj_ID filesep NF_dir filesep '_misc' filesep 'intersect_mask_NF_' subj_ID '_' tmp_normalize_prefix reg_vols_prefix(1:end-1) '.nii'],'file')
        fprintf(['\nRunning MNI intersect mask creation (for ' tmp_normalize_prefix reg_vols_prefix '*.nii files) for ' subj_ID '...\n']);
        intersect_mask = dir([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep tmp_normalize_prefix reg_vols_prefix '*000006-01.nii']); % again, do 6:265 to have 1:260 volumes out of 1:270
        if isempty(intersect_mask)
            error(['Unable to find 6th ' tmp_normalize_prefix reg_vols_prefix '*.nii volume in V01/Run1 folder for ' subj_ID '! Please double-check...']);
        else
            intersect_mask = load_untouch_nii([procpath subj_ID filesep NF_dir filesep 'V01' filesep 'Run1' filesep intersect_mask(1).name]);
            intersect_mask_data = double(intersect_mask.img);
            intersect_mask_data = intersect_mask_data(:);
            intersect_mask.fileprefix = ['intersect_mask_NF_' subj_ID '_' tmp_normalize_prefix reg_vols_prefix(1:end-1)];
        end
        for v = 1:length(all_NF_visits)
            all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
            for r = 1:length(all_NF_runs_per_visit)
                tmp_all_n_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name ...
                    filesep all_NF_runs_per_visit(r).name filesep tmp_normalize_prefix reg_vols_prefix '*.nii']);
                if ~isequal(length(tmp_all_n_files),260) && ~isequal(length(tmp_all_n_files),270)
                    error(['Incorrect number of ' tmp_normalize_prefix reg_vols_prefix '*.nii files found in folder ' ...
                        all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name ' for ' subj_ID '! Please check manually!']);
                end
                for f = 1:length(tmp_all_n_files)
                    tmp_vol_for_erosion = load_untouch_nii(...
                        [procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                        all_NF_runs_per_visit(r).name filesep tmp_all_n_files(f).name]);
                    tmp_vol_for_erosion = reshape(double(tmp_vol_for_erosion.img),numel(tmp_vol_for_erosion.img),1);
                    %intersect_mask_data = max(min(intersect_mask_data,tmp_s6uf_vol_for_erosion),0); % eroding the mask iteratively
                    intersect_mask_data = min(intersect_mask_data,tmp_vol_for_erosion); % [no need for max(.,0)]
                end
            end
        end
        intersect_mask.img = int16(reshape(intersect_mask_data,size(tmp_vol_for_erosion,1),size(tmp_vol_for_erosion,2),size(tmp_vol_for_erosion,3)));
        save_untouch_nii(intersect_mask,[procpath subj_ID filesep NF_dir filesep '_misc' filesep intersect_mask.fileprefix '.nii']);
        clear tmp_vol_for_erosion intersect_mask_data tmp_all_n_files all_NF_runs_per_visit;
        fprintf(['\nMNI intersect mask (' intersect_mask.fileprefix ') creation completed!\n']);
    else
        fprintf(['\nSkipping MNI intersect mask creation (for ' tmp_normalize_prefix reg_vols_prefix '*.nii volumes) for ' subj_ID ' (mask already exists!)...\n']);
    end
else
    fprintf(['\nSkipping MNI intersect mask creation for normalized ' reg_vols_prefix '*.nii volumes for ' subj_ID ' (nflag = 0).\n']); %#ok
end




%   11) Smoothing of wReg2_* (normalized, MNI) volumes;

if nflag
    fprintf(['\nRunning MNI smoothing (and masking) batch (NF) for ' subj_ID '...\n']);
    for v = 1:length(all_NF_visits)
        all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep 'Run*']);
        for r = 1:length(all_NF_runs_per_visit)
            tic;
            all_NF_Reg_files = dir([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                all_NF_runs_per_visit(r).name filesep tmp_normalize_prefix reg_vols_prefix '*.nii']);
            if ~isequal(length(all_NF_Reg_files),270) % && ~isequal(length(all_Reg_files),260)
                error(['Incorrect number of ' tmp_normalize_prefix reg_vols_prefix '*.nii files found in directory ' all_NF_visits(v).name ' ' ...
                    all_NF_runs_per_visit(r).name ' for ' subj_ID '! Found ' num2str(length(all_NF_Reg_files)) ' files instead of 270. Please check manually.']);
            end
            if exist([procpath subj_ID filesep NF_dir filesep all_NF_visits(v).name filesep ...
                    all_NF_runs_per_visit(r).name filesep 's' num2str(smoothing_factor) all_NF_Reg_files(1).name],'file')
                fprintf(['\nSkipping MNI smoothing for ' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ' (already exists)...\n']);
            else
                clear matlabbatch;
                for f = 1:length(all_NF_Reg_files)
                    matlabbatch{1}.spm.spatial.smooth.data{f,1} = [procpath subj_ID filesep NF_dir filesep ...
                        all_NF_visits(v).name filesep all_NF_runs_per_visit(r).name filesep all_NF_Reg_files(f).name ',1'];
                end
                matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = ['s' num2str(smoothing_factor)]; % [but no masking here]
                %spm('defaults','fmri');
                %spm_jobman('initcfg');
                fprintf(['\nSmoothing MNI ' tmp_normalize_prefix reg_vols_prefix '* volumes... (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ')\n']); % tic;
                spm_jobman('run',matlabbatch);
                fprintf(['Done (' subj_ID ' ' all_NF_visits(v).name ' ' all_NF_runs_per_visit(r).name ', ' num2str(toc) ' s).\n']);
            end
        end
    end
else
    fprintf(['\nSkipping MNI smoothing of ' reg_vols_prefix '*.nii volumes for ' subj_ID ' (nflag = 0).\n']); %#ok
end




% 
% %
% %%% #10 Smoothing of MNI-normalized volumes (write s6wReg_*.nii files)
% %
% %tmp_skip = true;
% 
% if ~tmp_skip
%     fprintf(['\nRunning MNI smoothing (and masking) batch (NF) for ' subj_ID '...\n']); %tic;
%     for v = 1:length(all_visits_dirs)
%         all_NF_runs_per_visit = dir([procpath subj_ID filesep NF_dir filesep all_visits_dirs(v).name filesep 'Run*']);
%         for r = 1:length(all_NF_runs_per_visit)
%             all_NF_Reg_files = dir([procpath subj_ID filesep NF_dir filesep all_visits_dirs(v).name filesep ...
%                 all_NF_runs_per_visit(r).name filesep 'wReg_*.nii']);
%             if ~isequal(length(all_NF_Reg_files),270) % && ~isequal(length(all_Reg_files),260)
%                 error(['Incorrect number of wReg_*.nii files found in directory ' all_visits_dirs(v).name ' ' ...
%                     all_NF_runs_per_visit(r).name ' for ' subj_ID '! Found ' num2str(length(all_NF_Reg_files)) ' files instead of 270. Please check manually.']);
%             end
%             if exist([procpath subj_ID filesep NF_dir filesep all_visits_dirs(v).name filesep ...
%                     all_NF_runs_per_visit(r).name filesep 's6' all_NF_Reg_files(1).name],'file')
%                 fprintf(['\nSkipping smoothing for ' subj_ID ' ' all_visits_dirs(v).name ' ' all_NF_runs_per_visit(r).name ' (already exists)...\n']);
%             else
%                 clear matlabbatch;
%                 for v = 1:length(all_NF_Reg_files)
%                     matlabbatch{1}.spm.spatial.smooth.data{v,1} = [procpath subj_ID filesep NF_dir filesep ...
%                         all_visits_dirs(v).name filesep all_NF_runs_per_visit(r).name filesep all_NF_Reg_files(v).name ',1'];
%                 end
%                 matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
%                 matlabbatch{1}.spm.spatial.smooth.dtype = 0;
%                 matlabbatch{1}.spm.spatial.smooth.im = 0;
%                 matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
%                 %spm('defaults','fmri');
%                 %spm_jobman('initcfg');
%                 %fprintf(['\nSmoothing... (' subj_ID ' ' all_visits_dirs(v).name ' ' all_runs_for_given_visit(r).name ')\n']); % tic;
%                 spm_jobman('run',matlabbatch);
%                 %fprintf(['Done (' subj_ID ' ' all_visits_dirs(v).name ' ' all_runs_for_given_visit(r).name ', ' num2str(toc) ' s)\n']);
%             end
%         end
%     end
% end




fprintf(['\nEVERYTHING DONE (' subj_ID ' NF pre-processing).\n']);




end



