%   Project Title: Mind In Motion - Age-differences uneven terrain walking
%
%   Code Designer: Jacob salminen & Chang Liu 
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

%   Chang Liu - 20231218
%   Epoch for resting/walking trials for all ICS
%   Done - Cleanup for publish - 20250930

%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
tic

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
%- define the directory to the src folderd
run_dir = [source_dir filesep 'scripts' filesep REPO_NAME filesep '2_GLOBAL_BATCH_CL'];
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
ADD_CLEANING_SUBMODS = 1;
setWorkspace_CL
%% PARPOOL SETUP ======================================================= %%
if ~ispc
    %## NOTE, you will need to edit icadefs's EEGOPTION_FILE to contain the
    %unix and pc paths for the option file on the M drive otherwise it just
    %does weird stuff. 
    pop_editoptions('option_storedisk', 1, 'option_savetwofiles', 1, ...
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
    pPool = parpool(pp, SLURM_POOL_SIZE, 'IdleTimeout', 1440);
else
    pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
    SLURM_POOL_SIZE = 1;
end
%% (PARAMETERS) ======================================================== %%
%## hard define
%- datset name
DATA_SET = 'STUDY-MiM_CRUNCH_202312';
%- study group and saving
SESSION_NUMBER = '1';
SAVE_ALLEEG = false;
SAVE_EEG = true; %true;
OVERRIDE_DIPFIT = true;
%- epoching params
DO_SLIDING_WINDOW = true;
%* sliding window
WINDOW_LENGTH = 5.25;   % sliding window length in seconds
PERCENT_OVERLAP = 0.0; % percent overlap between epochs % Chang Liu 12/18/2023 set to 0% to enough number of epochs
%* gait
EVENT_CHAR = 'RHS'; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
STD_TIMEWARP = 3;
EPOCH_TIME_LIMITS = [-1,4.25]; 
% with frequency decomposition artifact during ERSP creation
% paper
TIMEWARP_EVENTS = {'RHS', 'LTO', 'LHS', 'RTO', 'RHS'};
if DO_SLIDING_WINDOW
    SUFFIX_PATH_EPOCHED = 'SLIDING_EPOCHED_ALL';
%     TRIAL_TYPES = {'rest','0p25','0p5','0p75','1p0','flat','low','med','high'};
    TRIAL_TYPES = {'rest'};
else
    SUFFIX_PATH_EPOCHED = 'GAIT_EPOCHED_ALL';
    TRIAL_TYPES = {'0p25','0p5','0p75','1p0','flat','low','med','high'};
end
%- eeglab_cluster.m spectral params
% FREQ_LIMITS = [1,100];
% CYCLE_LIMITS = [3,0.8];
% SPEC_MODE = 'psd'; %'fft'; %'psd'; %options: 'psd','fft','pburg','pmtm'
% FREQ_FAC = 4;
% PAD_RATIO = 2;

%- datetime override
dt = '12012023_OAYA104_icc0p65-0p4_changparams';
%- Subject Directory information
OA_PREP_FPATH = '11262023_YAOAN104_iccRX0p65_iccREMG0p4_changparams'; 
%% DEFINE PATHS
%## soft define
STUDIES_DIR = [source_dir filesep DATA_SET filesep '_studies'];
OUTSIDE_DATA_DIR = [source_dir filesep DATA_SET filesep OA_PREP_FPATH]; % JACOB,SAL(02/23/2023)

DATA_DIR = RAW_DATA_DIR;

study_fName_1 = 'all_comps_study';
study_fName_2 = 'epoch_study_all';

save_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% Store fNames and fPaths
conditions      = cell(1,length([SUBJ_ITERS{:}]));
groups          = cell(1,length([SUBJ_ITERS{:}]));
sessions        = cell(1,length([SUBJ_ITERS{:}]));
subjectNames    = cell(1,length([SUBJ_ITERS{:}]));
fNames          = cell(1,length([SUBJ_ITERS{:}]));
fPaths          = cell(1,length([SUBJ_ITERS{:}]));
chanlocs_fPaths   = cell(1,length([SUBJ_ITERS{:}]));
dipfit_norm_fPaths = cell(1,length([SUBJ_ITERS{:}]));
stack_iter = 0;
for group_i = 1:length(SUBJ_ITERS)
    sub_idx = SUBJ_ITERS{group_i}; %1:2; %1:length(SUBJ_PICS{GROUP_INT}); %1:2;
    %- set cnt
    cnt = stack_iter + 1;
    %## Assigning paths for .set, headmodel,& channel file
    for subj_i = sub_idx
        %- ICA fPaths
        fPaths{cnt} = [OUTSIDE_DATA_DIR filesep SUBJ_PICS{group_i}{subj_i} filesep 'clean'];
%         fPaths{cnt} = [load_dir filesep SUBJ_PICS{group_i}{subj_i} filesep 'ICA'];
        tmp = dir([fPaths{cnt} filesep '*.set']);
        try
            fNames{cnt} = tmp.name;
            %- Chanlocs fPaths
            chanlocs_fPaths{cnt} = [DATA_DIR filesep DATA_SET filesep SUBJ_PICS{group_i}{subj_i} filesep 'MRI' filesep 'CustomElectrodeLocations.mat'];
            dipfit_norm_fPaths{cnt} = [fPaths{cnt} filesep 'dipfit_fem_norm_ants.mat'];% use normalized dipfit folder
            %- Prints
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('ICA Exists: %i\n',(exist([fPaths{cnt} filesep fNames{cnt}],'file') && exist([fPaths{cnt} filesep 'W'],'file')))
            fprintf('Normalized DIPFIT Exists: %i\n',exist(dipfit_norm_fPaths{cnt},'file'));
        catch e
            fprintf('==== Subject %s Paths ====\n',SUBJ_PICS{group_i}{subj_i})
            fprintf('%s\n',getReport(e))
            dipfit_norm_fPaths{cnt} = [];
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
inds = logical(cellfun(@(x) exist(x,'file'),dipfit_norm_fPaths));
chanlocs_fPaths = chanlocs_fPaths(inds);
dipfit_norm_fPaths = dipfit_norm_fPaths(inds);
fPaths = fPaths(inds);
fNames = fNames(inds);
sessions = sessions(inds);
groups = groups(inds);
conditions = conditions(inds);
subjectNames = subjectNames(inds);

%% INITIALIZE PARFOR LOOP VARS
if exist('SLURM_POOL_SIZE','var')
%     POOL_SIZE = min([SLURM_POOL_SIZE,length(MAIN_ALLEEG)]);
    POOL_SIZE = min([SLURM_POOL_SIZE,length(fPaths)]);
else
    POOL_SIZE = 1;
end
% fPaths = {MAIN_ALLEEG.filepath};
% fNames = {MAIN_ALLEEG.filename};
LOOP_VAR = 1:length(fPaths);
%- clear vars for memory
% clear MAIN_ALLEEG
%% GENERATE EPOCH MAIN FUNC
%## PARFOR LOOP
parfor (subj_i = LOOP_VAR,POOL_SIZE)
    %## LOAD EEG DATA
    EEG = mim_create_alleeg(fNames(subj_i),fPaths(subj_i),subjectNames(subj_i),save_dir,...
                        conditions(subj_i),groups(subj_i),sessions(subj_i));
%     EEG = MAIN_ALLEEG(subj_i);
%     EEG = pop_loadset('filepath',fPaths{subj_i},'filename',fNames{subj_i});
    
    %%-------------- Reject eye component
    fprintf('Rejecting eye movement...\n');
    EEG = reject_eye_ic(EEG);

    fprintf('Running subject %s\n',EEG.subject)
    %- Recalculate ICA Matrices && Book Keeping
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end
    
    %## PARSE TRIALS
    epoched_fPath = [save_dir filesep EEG.subject filesep SUFFIX_PATH_EPOCHED];
    fPath = [epoched_fPath filesep [TRIAL_TYPES{:}]];
    fName = sprintf('%s_%s_EPOCH_TMPEEG_ALL.set',EEG.subject,[TRIAL_TYPES{:}]);
    if ~exist(fPath,'dir')
        mkdir(fPath)
    end
    %- parse
    try
        %## EPOCH
        [ALLEEG,timewarp_struct] = mim_parse_trials(EEG,DO_SLIDING_WINDOW,...
            'EPOCH_TIME_LIMITS',EPOCH_TIME_LIMITS,...
            'STD_TIMEWARP',STD_TIMEWARP,...
            'COND_CHARS',TRIAL_TYPES,...
            'WINDOW_LENGTH',WINDOW_LENGTH,...
            'PERCENT_OVERLAP',PERCENT_OVERLAP);%add by CL
        %## REMOVE USELESS EVENT FIELDS (Improve Load Time)
        for i = 1:length(ALLEEG)
            if isfield(ALLEEG(i).event,'trialName')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'trialName');
            end
            if isfield(ALLEEG(i).event,'channel')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'channel');
            end
            if isfield(ALLEEG(i).event,'code')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'code');
            end
            if isfield(ALLEEG(i).event,'bvtime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvtime');
            end
            if isfield(ALLEEG(i).event,'bvmknum')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'bvmknum');
            end
            if isfield(ALLEEG(i).event,'datetime')
                ALLEEG(i).event = rmfield(ALLEEG(i).event,'datetime');
            end
        end
        %## SAVE EEG's AS INDIVIDUAL FILES (CONNECTIVITY)
        cond_files = struct('fPath',[],'fName',[]);
        if SAVE_ALLEEG
            for i = 1:length(ALLEEG)
                %- save each parsed trial/condition to own folder to help save
                %memory. EEGLAB is weird like that.
                REGEX_FNAME = 'cond_%s';
                tmp_fPath = [epoched_fPath filesep sprintf(REGEX_FNAME,ALLEEG(i).condition)];
                if ~exist(tmp_fPath,'dir')
                    mkdir(tmp_fPath)
                end
                [~] = pop_saveset(ALLEEG(i),'savemode','twofiles',...
                    'filepath',tmp_fPath,'filename',sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition));
                cond_files(i).fPath = tmp_fPath;
                cond_files(i).fName = sprintf([REGEX_FNAME '.set'],ALLEEG(i).condition);
            end
        end
        ALLEEG = pop_mergeset(ALLEEG,1:length(ALLEEG),1);
        ALLEEG.etc.cond_files = cond_files;
        %## timewarp for across condition
        if ~DO_SLIDING_WINDOW
            timewarp = make_timewarp(ALLEEG,TIMEWARP_EVENTS,'baselineLatency',0, ...
                    'maxSTDForAbsolute',inf,...
                    'maxSTDForRelative',inf);
            %- subject specific warpto (later use to help calc grand avg warpto across subjects)
            timewarp.warpto = nanmedian(timewarp.latencies);        
            goodepochs  = sort([timewarp.epochs]);
            %- probably not needed? 
            sedi = setdiff(1:length(ALLEEG.epoch),goodepochs);
            %- reject outlier strides
            ALLEEG = pop_select(ALLEEG,'notrial',sedi);
            %- store timewarp structure in EEG
            ALLEEG.timewarp = timewarp;
    %         disp(EEG.subject); disp(allWarpTo); disp(grandAvgWarpTo);
            %- store condition-by-conditino timewarpings
            ALLEEG.etc.timewarp_by_cond = timewarp_struct;
            %## STRUCT EDITS
            ALLEEG.urevent = []; % might be needed
            ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
        end
        %## STRUCT EDITS
        ALLEEG.urevent = []; % might be needed
        ALLEEG.etc.epoch.epoch_limits = EPOCH_TIME_LIMITS;
        %- checks
        ALLEEG = eeg_checkset(ALLEEG,'eventconsistency');
        ALLEEG = eeg_checkset(ALLEEG);
        ALLEEG = eeg_checkamica(ALLEEG);
        %- save
        pop_saveset(ALLEEG,'savemode','twofiles',...
                'filename',fName,...
                'filepath',fPath,...
                'version','6');
%         tmp{subj_i} = ALLEEG;
    catch e
        rmv_subj(subj_i) = 1;
        EEG.timewarp = struct([]);
        EEG.urevent = [];
%         tmp{subj_i} = []; %EEG;
        fprintf(['error. identifier: %s\n',...
                 'error. %s\n',...
                 'error. on subject %s\n',...
                 'stack. %s\n'],e.identifier,e.message,EEG.subject,getReport(e));
    end
    
    %% To save memory, clean unnecessary files
    EEG = []; ALLEEG = []; 
end
