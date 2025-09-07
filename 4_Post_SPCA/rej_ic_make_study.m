% Create study for cluster ICs after rejecting ICs that are not brain 

% Chang Liu - 2021-11-23 - V1
% Chang Liu - 2025-5-7 Find a bug that the code is adding wrong comps to
% the subject. Add code to check if the subject name match

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.


%   NJacobsen notes
%   When timewarping data, save values as EEG.timewarp = timewarp;
%   EEG.timewarp.medianlatency = median(timewarp.latencies(:,:));%Warping to the median latency of my 5 events
%   By default, std_ersp will use the median of all subject's
%   timewarp.latencies(:,:) as 'timewarpms' unless individual subject 
%   warpto is indiciated using 'timewarpms', 'subject tw matrix'
%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230417.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

%{
%## RESTORE MATLABs
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'liu.chang1'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);
%- CD
% cfname_path    = mfilename('fullpath');
% cfpath = strsplit(cfname_path,filesep);
% cd(cfpath);
%% (EDIT: PATH TO YOUR GITHUB REPO) ==================================== %%
%- GLOBAL VARS
REPO_NAME = 'MiM_CRUNCH';
%- determine OS
if strncmp(computer,'PC',2)
    PATH_ROOT = ['M:' filesep USER_NAME filesep 'scripts' filesep REPO_NAME]; % path 2 your github folder
    source_dir = ['M:' filesep USER_NAME ];
    
    %########### USE JACOB's Folder for RAW_DATA
    RAW_DATA_DIR = ['M:' ...
        filesep 'jsalminen' filesep 'GitHub' filesep 'par_EEGProcessing' filesep 'src' ...
        filesep '_data' filesep 'MiM_dataset'];

    %##############################################

else  % isunix
    PATH_ROOT = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME filesep 'scripts' filesep REPO_NAME]; % path 2 your github folder
    source_dir = [filesep 'blue' filesep 'dferris',...
        filesep USER_NAME ];
    
    %########### USE JACOB's Folder for RAW_DATA
    RAW_DATA_DIR = [filesep 'blue' filesep 'dferris',...
        filesep 'jsalminen' filesep 'GitHub' filesep 'par_EEGProcessing' filesep 'src' ...
        filesep '_data' filesep 'MiM_dataset'];

    %##############################################

end
run_dir = [source_dir filesep 'scripts' filesep REPO_NAME filesep '4_Post_SPCA'];
%% CD ================================================================== %%
%- cd to run directory
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
addpath(PATH_ROOT)

% -- Run config file
MiM_CRUNCH_config
%% SET WORKSPACE ======================================================= %%
global ADD_CLEANING_SUBMODS
ADD_CLEANING_SUBMODS = false;
setWorkspace_CL
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0,'option_saveversion6',1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    disp(['SLURM_JOB_ID: ', getenv('SLURM_JOB_ID')]);
    disp(['SLURM_CPUS_ON_NODE: ', getenv('SLURM_CPUS_ON_NODE')]);
    %## allocate slurm resources to parpool in matlab
    %- get cpu's on node and remove a few for parent script.
    SLURM_POOL_SIZE = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %- create cluster
    pp = parcluster('local');
    %- Number of workers for processing (NOTE: this number should be higher
    %then the number of iterations in your for loop)
    fprintf('Number of workers: %i\n',pp.NumWorkers);
    fprintf('Number of threads: %i\n',pp.NumThreads);
    %- make meta data dire1ory for slurm
    mkdir([run_dir filesep getenv('SLURM_JOB_ID')])
    pp.JobStorageLocation = strcat([run_dir filesep], getenv('SLURM_JOB_ID'));
    %- create your p-pool (NOTE: gross!!)
    pPool = parpool(pp, SLURM_POOL_SIZE);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (JS_PARAMETERS) ===================================================== %%
%## hard define
%- datset name
DATA_SET = 'STUDY-MiM_CRUNCH_202312';
TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high','rest'};
if ispc
    FIELDTRIP_PATH = 'M:\jsalminen\GitHub\par_EEGProcessing\submodules\fieldtrip';
else
    FIELDTRIP_PATH = '/blue/dferris/jsalminen/GitHub/par_EEGProcessing/submodules/fieldtrip';
end
%- compute measures for spectrum and ersp
FORCE_RECALC_SPEC = true;
%- statistics & conditions
speed_trials = {'0p25','0p5','0p75','1p0','rest'};
terrain_trials = {'flat','low','med','high','rest'};
COND_DESIGNS = {speed_trials,terrain_trials};
SESSION_NUMBER = '1';
%##
%* ERSP PARAMS
ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/16/2023) JS, updating method to bootstrap as per CL YA paper
% (07/19/2023) JS, subbaseline set to off 
% (07/20/2023) JS, subbaseline set to on, generates different result? 
SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);
% (07/16/2023) JS, updating mcorrect to fdr as per CL YA paper
% (07/31/2023) JS, changing fieldtripnaccu from 2000 to 10000 to match CL's
% pipeline although this doesn't align with her YA manuscript methods?
SPEC_PARAMS = struct('freqrange',[1,200],...
    'subject','',...
    'specmode','psd',...
    'freqfac',4,...
    'logtrials','on',...
    'comps','all',...
    'plot_freqrange',[4,60],...
    'plot_ylim',[-35,-8],...
    'subtractsubjectmean','on',...
    'plotmode','normal');
ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
% NOTE: see. statcondfieldtrip.m or std_stat.m
COND_EVENT_CHAR = 'cond';
%- clustering parameters
MIN_ICS_SUBJ = [5]; %[2,3,4,5,6,7,8]; % iterative clustering
% K_RANGE = [10,22];
MAX_REPEATED_ITERATIONS = 1;
CLUSTER_SWEEP_VALS = [12,14,19]; %[10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
%***********************************************
% DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = true;
% DO_K_SWEEPING = false;
%************************************************
% (08/21/2023) JS, this currenlty doesn't do anything but take up more
% memory.
REPEATED_CLUSTERING_STD = 3;
CLUSTER_PARAMS = struct('algorithm','kmeans',...
    'clust_num',20,...
    'save','off',...
    'filename',[],...
    'filepath',[],...
    'outliers',inf(),...
    'maxiter',200);
%- 
%- custom params
% colormap_ersp = othercolor('RdYlBu11');
% colormap_ersp = colormap_ersp(end:-1:1,:);
%NOTE: (NJacobsen); warp each subject's tw matrix to the entire group's median event
%latencies [1=ON], or use individual subject's median event latencies [0=OFF]. TW must be ON
%for this setting to do anything.
clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = 'dipole_1'; %['dipole_',num2str(clustering_weights.dipoles),...
%     '_scalp_',num2str(clustering_weights.scalp),'_ersp_',num2str(clustering_weights.ersp),...
%     '_spec_',num2str(clustering_weights.spec)];
STD_PRECLUST_COMMAND = {'dipoles','weight',clustering_weights.dipoles};
%- iterative clustering parameters
% n_iterations = 50;
% outlier_sigma = 3;
%- datetime override
%- datetime override
dt = '12012023_OAYA104_icc0p65-0p4_changparams';
%- Subject Directory information
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams'; 
%## soft define
DATA_DIR = [source_dir ];
OUTSIDE_DATA_DIR = [DATA_DIR filesep DATA_SET filesep OA_PREP_FPATH];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);

study_fName_1 = 'all_comps_study';
study_fName_2 = 'epoch_study';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];

save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% Store fNames and fPaths
conditions      = cell(1,length([SUBJ_ITERS{:}]));
groups          = cell(1,length([SUBJ_ITERS{:}]));
sessions        = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames_epoch          = cell(1,length([SUBJ_ITERS{:}]));
fPaths_epoch          = cell(1,length([SUBJ_ITERS{:}]));
chanlocs_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
% dipfit_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = cell(1,length([SUBJ_ITERS{:}]));
% vol_fPaths = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        path_epoch = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'GAIT_EPOCHED_ALL' filesep '0p250p50p751p0flatlowmedhigh'];
        if ~isempty(dir([path_epoch filesep '*_addinfo.set']))
%             tmp = dir([fPaths{cnt} filesep '*.set']);
            tmp_epoch = dir([path_epoch filesep '*.set']);
            try
%                 fNames{cnt} = tmp.name;
                fNames_epoch{cnt} = tmp_epoch.name;
                fPaths_epoch{cnt} = path_epoch;

                %- Chanlocs fPaths
        %         chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'EEG' filesep 'HeadScan' filesep 'CustomElectrodeLocations.mat'];
%                 chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
%         %         dipfit_fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'head_model' filesep 'dipfit_struct.mat'];
%     %             dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm.mat'];
%                 dipfit_norm_fPaths{cnt} = [fPaths_epoch{cnt} filesep 'dipfit_fem_norm_ants.mat'];
%                 %- Prints
                fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
%                 fprintf('ICA Exists: %i\n',(exist([fPaths_epoch{cnt} filesep fNames{cnt}],'file') && exist([fPaths_epoch{cnt} filesep 'W'],'file')))
%         %         fprintf('DIPFIT Exists: %i\n',exist(dipfit_fPaths{cnt},'file'));
%                 fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{cnt},'file'));
            catch e
                fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
                fprintf('%s\n',getReport(e))
                
            end
        else 
            fprintf('No file found for %s \n',SUBJ_PICS{group_i}{subj_i})
        end
        cnt = cnt + 1;
    end
    %- reset cnt
    cnt = stack_iter + 1;
    %## Assigning paths for eeglab study
    for subj_i = sub_idx
        subjectNames{cnt} = SUBJ_PICS{group_i}{subj_i};
        tmp = join(TRIAL_TYPES,'_'); 
        conditions{cnt} = tmp{:};
        groups{cnt} = GROUP_NAMES{group_i};
        sessions{cnt} = SESSION_NUMBER;
        cnt = cnt + 1;
    end
    stack_iter = stack_iter + length(SUBJ_ITERS{group_i});
end
%- remove subjects without a dipole fit
inds = logical(cellfun(@(x) exist(x,'file'),fPaths_epoch));
% chanlocs_fPaths = chanlocs_fPaths(inds);
% dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
fPaths_epoch = fPaths_epoch(inds);
fNames_epoch = fNames_epoch(inds);
% sessions = sessions(inds);
% groups = groups(inds);
% conditions = conditions(inds);
subjectNames = subjectNames(inds);

%% 
LOOP_VAR = 1:length(fPaths_epoch);

%% Make study and reject ICs
OVERWRITE_study = true;
parfor subj_i = LOOP_VAR
    save_dir = fullfile(fPaths_epoch{subj_i},'reject_ic');
    if ~exist(fullfile(save_dir,[subjectNames{subj_i},'_ICRej.mat']),'file')|OVERWRITE_study
        EEG = pop_loadset('filename',fNames_epoch{subj_i},'filepath',fPaths_epoch{subj_i});
        if ~exist(save_dir,'dir');mkdir(save_dir);end
        [EEG] = mim_reject_ics_CL(EEG,save_dir);
    else
        fprintf('Skip %s \n',subjectNames{subj_i});
        continue;
    end
end                

