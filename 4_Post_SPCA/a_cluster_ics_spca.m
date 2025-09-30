% Create study for cluster ICs. This code only works for cluster without
% using ERSP. Precompute ERSP needed to be done on Hipergator
% Chang Liu - 2021-11-23 - V1
% Chang Liu - 2025-05-06 - Find an error in clustering outcome, some
% participants' comps do not match with their rejection IC.mat. Really
% weird mistake.

% Chang Liu - clean up for the published code 

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%   Code Designer: Jacob salminen, Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20230417.0
%   Previous Version: n/a
%   Summary: The following script is to identify potential brain components
%   for the Mind-In-Motion study

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_a_cluster_ics.sh


%% (REQUIRED SETUP 4 ALL SCRIPTS) ====================================== %%
%- DATE TIME
dt = datetime;
dt.Format = 'MMddyyyy';
%- VARS
USER_NAME = 'liu.chang1'; %getenv('username');
fprintf(1,'Current User: %s\n',USER_NAME);

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
run_dir = [source_dir filesep 'scripts' filesep REPO_NAME filesep '3_ANALYZE_CL'];
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

SPEC_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','perm',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','fdr',...  % ['cluster'|'fdr']
    'fieldtripnaccu',2000);

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
MIN_ICS_SUBJ = [5]; % Minimal brain IC for each participant iterative clustering 

MAX_REPEATED_ITERATIONS = 1;
CLUSTER_SWEEP_VALS = [10,11]; %k = [9:14] I sweep from 9 - 14 clusters
reduce_method = 'max_icalabel'; %The way to reduce to 1 IC per person per cluster
%***********************************************
% DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = 1;% To remove people with fewer than 5 brain components 
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

clustering_weights.dipoles = 1;
clustering_weights.scalp = 0;
clustering_weights.ersp = 0;
clustering_weights.spec = 0;
clustering_method = 'dipole_1'; 
STD_PRECLUST_COMMAND = {'dipoles','weight',clustering_weights.dipoles};
%- iterative clustering parameters
% n_iterations = 50;
% outlier_sigma = 3;
%- datetime override
dt = '12012023_OAYA104_icc0p65-0p4_changparams';
%## soft define
DATA_DIR = [source_dir ];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);

study_fName_1 = 'rest_study_new_spca_all_fixed_reduce_202505';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];

save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%% Make the giant study
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
OVERWRITE_study = true;
if ~exist([load_dir filesep study_fName_1 '.study'],'file')|OVERWRITE_study    
    % Create study - Chang Liu
    [STUDY,ALLEEG] =  mim_create_epoch_study_rej_ic([],'epoch_study_spca_all_fixed_202505',load_dir,...
    ['GAIT_EPOCHED_ALL' filesep '0p250p50p751p0flatlowmedhigh'],'GAIT_EPOCHED_ALL','SLIDING_EPOCHED_ALL');
    [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
                                            STUDY.filename,STUDY.filepath,...
                                            'RESAVE_DATASETS','off');                        
    [STUDY_rest,ALLEEG_rest] = mim_create_epoch_study_rej_ic([],study_fName_1,load_dir,['SLIDING_EPOCHED_ALL' filesep 'rest'],...
        'GAIT_EPOCHED_ALL','SLIDING_EPOCHED_ALL');   
    [STUDY_rest,ALLEEG_rest] = parfunc_save_study(STUDY_rest,ALLEEG_rest,...
                                            STUDY_rest.filename,STUDY_rest.filepath,...
                                            'RESAVE_DATASETS','off');  
%     if DO_K_ICPRUNE                   
%         STUDY = STUDY_rest; ALLEEG = ALLEEG_rest;
%     end
else
    if ~ispc
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    else
        [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
    end
end

CLUSTER_PARAMS.filename = STUDY.filename;
CLUSTER_PARAMS.filepath = STUDY.filepath;


%% SETUP the Timewarp parameters
% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
%     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
% grandAvgWarpTo = floor(nanmean(allWarpTo)); % tends to be longer? (e.g., [0,262,706,982,1415])
grandAvgWarpTo = [0         256         728         999        1461]; % set to be number to ease computation

% ---------------------------------------------------------------------------
% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
% b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)-1000*(1/ALLEEG(1).srate)];
b_lims =[grandAvgWarpTo(1) grandAvgWarpTo(5)];
% ERSP_CROP_TIMES=[grandAvgWarpTo(1)+abs(ALLEEG(1).etc.epoch.epoch_limits(1))*1000, grandAvgWarpTo(5)];
ERSP_CROP_TIMES=[grandAvgWarpTo(1), grandAvgWarpTo(5)];
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',b_lims(1),b_lims(2));
disp(grandAvgWarpTo);
%% (SET PARAMS)
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
tmp_group_orig = cell(length(ALLEEG),1);
tmp_group_unif = cell(length(ALLEEG),1);
for subj_i = 1:length(ALLEEG)
    tmp_group_orig{subj_i} = ALLEEG(subj_i).group;
    tmp_group_unif{subj_i} = 'Older Adults';
end
%-
%{
for subj_i = 1:length(ALLEEG)
    ALLEEG(subj_i).group = tmp_group_orig{subj_i};
    STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
end
%}
%- NOTE: partly adapt from bemobil_repeated_clustering
numIC = zeros(length(STUDY.datasetinfo),1);
for n = 1:length(STUDY.datasetinfo)
    numIC(n) = size(STUDY.datasetinfo(n).comps,2);
end
fprintf('Mean Clusters = %s \n',num2str(mean(numIC)));
mean_IC_allSub = floor(mean(numIC)+10);
%% (PRECOMPUTE MEASURES) COMPUTE SPECTRUMS
DO_precompute = false; 
overwrite = 'on';
DO_TIMEWARP = true;

if DO_precompute
    tmp = strsplit(ALLEEG(1).filename,'.');
    spec_f = [ALLEEG(1).filepath filesep sprintf('%s.icaspec',ALLEEG(1).subject)];
    topo_f = [ALLEEG(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
    % pPool = parpool(pp, SLURM_POOL_SIZE);
    if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
        fprintf('Calculating Spectograms...\n');

        %- override variables for the stats
        for subj_i = 1:length(ALLEEG)
            ALLEEG(subj_i).group = tmp_group_orig{subj_i};
            STUDY.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
            for in_i = 1:length(STUDY.datasetinfo(subj_i).trialinfo)
                STUDY.datasetinfo(subj_i).trialinfo(in_i).group = tmp_group_orig{subj_i};
            end
        end
        for subj_i = 1:length(ALLEEG),ceil(length(ALLEEG)/2)
            EEG = ALLEEG(subj_i);
            TMP_STUDY = STUDY;
            EEG = eeg_checkset(EEG,'loaddata');
            if isempty(EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
            end
            %- overrride datasetinfo to trick std_precomp to run.
            TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
            fprintf('SUBJECT: %s\n',TMP_STUDY.datasetinfo.subject);
            fprintf('GROUP: %s\n',TMP_STUDY.datasetinfo.group);
            disp(STUDY.datasetinfo(subj_i));
            TMP_STUDY.datasetinfo(1).index = 1;
            [~, ~] = std_precomp(TMP_STUDY, EEG,...
                        'components',...
                        'recompute',overwrite,...
                        'spec','on',...
                        'scalp','on',...
                        'savetrials','on',...
                        'specparams',...
                        {'specmode',SPEC_PARAMS.specmode,'freqfac',SPEC_PARAMS.freqfac,...
                        'freqrange',SPEC_PARAMS.freqrange,'logtrials',SPEC_PARAMS.logtrials});
                                 
             fprintf('Computing ERSP...');
             
             if DO_TIMEWARP
                timewarp_param = EEG.timewarp.latencies;
                timewarpms_param = grandAvgWarpTo;
             else
                timewarp_param = [];
                timewarpms_param = [];
             end
             %-
             filepath = TMP_STUDY.datasetinfo(1).filepath;
             trialinfo = std_combtrialinfo(TMP_STUDY.datasetinfo, 1);
             filebase = fullfile(filepath, EEG.subject);
             
            [~, ~] = std_precomp(TMP_STUDY, EEG, 'comps', 'savetrials','on','recompute','on','ersp','on',...
            'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,'ntimesout',TIMEWARP_NTIMES,...
            'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'baseline',nan,...
            'timewarp',timewarp_param,'timewarpms',timewarpms_param},'itc','off'); %ERSP
        end
    end
    %--- Do the same for the rest condition
    tmp = strsplit(ALLEEG_rest(1).filename,'.');
    spec_f = [ALLEEG_rest(1).filepath filesep sprintf('%s.icaspec',ALLEEG_rest(1).subject)];
    topo_f = [ALLEEG_rest(1).filepath filesep sprintf('%s.icatopo',tmp{1})];
    % pPool = parpool(pp, SLURM_POOL_SIZE);
    if ~exist(spec_f,'file') || ~exist(topo_f,'file') || FORCE_RECALC_SPEC
        fprintf('Calculating Spectograms...\n');

        %- override variables for the stats
        for subj_i = 1:length(ALLEEG_rest)
            ALLEEG_rest(subj_i).group = tmp_group_orig{subj_i};
            STUDY_rest.datasetinfo(subj_i).group = tmp_group_orig{subj_i};
            for in_i = 1:length(STUDY_rest.datasetinfo(subj_i).trialinfo)
                STUDY_rest.datasetinfo(subj_i).trialinfo(in_i).group = tmp_group_orig{subj_i};
            end
        end
        for subj_i = 1:length(ALLEEG_rest),ceil(length(ALLEEG_rest)/2)
            EEG = ALLEEG_rest(subj_i);
            TMP_STUDY = STUDY_rest;
            EEG = eeg_checkset(EEG,'loaddata');
            if isempty(EEG.icaact)
                fprintf('%s) Recalculating ICA activations\n',EEG.subject);
                EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
                EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
            end
            %- overrride datasetinfo to trick std_precomp to run.
            TMP_STUDY.datasetinfo = STUDY_rest.datasetinfo(subj_i);
            fprintf('SUBJECT: %s\n',TMP_STUDY.datasetinfo.subject);
            fprintf('GROUP: %s\n',TMP_STUDY.datasetinfo.group);
            disp(STUDY_rest.datasetinfo(subj_i));
            TMP_STUDY.datasetinfo(1).index = 1;
            [~, ~] = std_precomp(TMP_STUDY, EEG,...
                        'components',...
                        'recompute',overwrite,...
                        'spec','on',...
                        'scalp','on',...
                        'savetrials','on',...
                        'specparams',...
                        {'specmode',SPEC_PARAMS.specmode,'freqfac',SPEC_PARAMS.freqfac,...
                        'freqrange',SPEC_PARAMS.freqrange,'logtrials',SPEC_PARAMS.logtrials});
        end
    end  
end
%% (STEP 4) CALCULATE CLUSTER SOLUTIONS AFTER PRUNING SUBJECTS WITH LOW NUMBER OF ICs
%## GET EEGLAB PATH
tmp = strsplit(path,';');
% tmp = strsplit(path,';');
b1 = regexp(tmp,'eeglab','end');
b2 = tmp(~cellfun(@isempty,b1));
PATH_EEGLAB = b2{1}; %(1:b1{1});
fprintf('EEGLAB path: %s\n',PATH_EEGLAB);
%- set default paths for boundary element head model
PATH_EEGLAB_BEM  = [PATH_EEGLAB filesep 'plugins' filesep 'dipfit' filesep 'standard_BEM' filesep];
% MNI_VOL = [PATH_EEGLAB_BEM filesep 'standard_vol.mat'];
%-
STUDY_COND_DESI = {{'subjselect',{},...
            'variable1','cond','values1',terrain_trials,...
            'variable2','group','values2',{}},...
            {'subjselect',{},...
            'variable1','cond','values1',speed_trials,...
            'variable2','group','values2',{}}};
        
HIRES_TEMPLATE = 'M:\liu.chang1\scripts\MiM_CRUNCH\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
if ~ispc
    HIRES_TEMPLATE = convertPath2UNIX(HIRES_TEMPLATE);
else
    HIRES_TEMPLATE = convertPath2Drive(HIRES_TEMPLATE);
end
%- assign hires_template default
tmp = strsplit(HIRES_TEMPLATE,filesep);
fpath = strjoin(tmp(1:end-1),filesep);
fname = tmp{end};
ext = strsplit(fname,'.');
fname = ext{1};
ext = ext{end};
MNI_VOL = [fpath filesep fname '_dipplotvol.mat'];
MNI_MRI = [fpath filesep fname '_dipplotmri.mat'];
% MNI_MESH = [fpath filesep fname '_mesh.mat'];
mri = load(MNI_MRI);
mri = mri.mri;
vol = MNI_VOL;
%% remove participants with low number of ICs and remove people with missing conditions
if DO_K_ICPRUNE 
    fprintf('performing pruning...\n');
%     for i = 1:length(MIN_ICS_SUBJ)
    for i = length(MIN_ICS_SUBJ)
        tmp_dir = [save_dir filesep sprintf('icrej_spca_%s_%i',study_fName_1,MIN_ICS_SUBJ(i))];
        if ~exist(tmp_dir,'dir')
            mkdir(tmp_dir)
        end
        ics_subjs = cellfun(@(x) size(x,2),{STUDY.datasetinfo.comps});
        idx1 = ics_subjs > MIN_ICS_SUBJ(i);%MIN_ICS_SUBJ decides the number of ICs as cutoff
        fprintf('removing %s participants \n',num2str(length(STUDY.datasetinfo)-sum(idx1)));
        
        disp(STUDY.subject(ics_subjs < MIN_ICS_SUBJ(i)));
%         TMP_ALLEEG = ALLEEG_rest(idx);
        
        condition_store = cell(length(STUDY.datasetinfo),8);%NH3007 doesn't complete high condition
        for ii = 1:length(STUDY.datasetinfo)
            condition_store(ii,1:length(unique({STUDY.datasetinfo(ii).trialinfo.cond}))) = unique({STUDY.datasetinfo(ii).trialinfo.cond});            
        end
        condition_store(cellfun(@isempty,condition_store)) = {'x'};
        idx2 = sum(contains(condition_store,{'flat','low','med','high'}),2) == 4;
        TMP_ALLEEG =  ALLEEG_rest(idx1 & idx2');    

        %##
        [TMP_STUDY, TMP_ALLEEG] = std_editset([],TMP_ALLEEG,...
                                        'updatedat','off',...
                                        'savedat','off',...
                                        'name',sprintf('temp_study_rejics_spca_%s_%i',study_fName_1,MIN_ICS_SUBJ(i)),...
                                        'filename',sprintf('temp_study_rejics_spca_%s_%i.study',study_fName_1,MIN_ICS_SUBJ(i)),...
                                        'filepath',tmp_dir);
                 
        % Copy the comps from the original STUDY to the TMP STUDY
        TMP_STUDY = rmfield(TMP_STUDY,'cluster');
        for j = 1:length(TMP_STUDY.datasetinfo)
            j_ind = find(strcmp({STUDY_rest.datasetinfo(:).subject},TMP_STUDY.datasetinfo(j).subject));
            TMP_STUDY.datasetinfo(j).comps = STUDY_rest.datasetinfo(j_ind).comps;
        end
%         [TMP_STUDY, TMP_ALLEEG] = std_editset(TMP_STUDY,TMP_ALLEEG,...
%                                        'commands',{'inbrain','on'});  
                                   
        [TMP_STUDY,TMP_ALLEEG] = std_checkset(TMP_STUDY,TMP_ALLEEG);
        
        %##
        [TMP_STUDY,TMP_ALLEEG] = std_preclust(TMP_STUDY,TMP_ALLEEG,1,STD_PRECLUST_COMMAND);
        
%         [TMP_STUDY,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY,TMP_ALLEEG,...
%                                             TMP_STUDY.filename,TMP_STUDY.filepath,...
%                                             'RESAVE_DATASETS','off');
        %- remove components outside the brain
%                        
        %- store essential info in STUDY struct for later reading
        freqrange = [4 60];
        TMP_STUDY.etc.clustering.preclustparams.clustering_weights = clustering_weights;
        TMP_STUDY.etc.clustering.preclustparams.freqrange = freqrange;
        all_solutions = mim_cluster_process(TMP_STUDY,TMP_ALLEEG,tmp_dir,...
            'CLUSTER_PARAMS',CLUSTER_PARAMS,...
            'REPEATED_CLUSTERING_STD',REPEATED_CLUSTERING_STD,...
            'MAX_REPEATED_ITERATIONS',MAX_REPEATED_ITERATIONS,...
            'KS_TO_TEST',CLUSTER_SWEEP_VALS);
        
        
        save_clustering_solutions = true;
        if save_clustering_solutions
            disp('Save and reduce clustering solutions');
        %       temporary
        %## (Step 2) CALCULATE REPEATED CLUSTERED SOLUTIONS
        %- NOTE: the clustering solutions are not exactly the same but not super different
    %     cluster_ks = K_RANGE(1):K_RANGE(2);
    %     parfor (j = 1:length(cluster_ks),length(cluster_ks))
        
            %}
            %% ## REMOVE DIPOLES OUTSIDE THE BRAIN AREA

            disp('RUNNING CLUSTER after removing dipoles outside brain');
            TMP_STUDY.etc.cluster_vars = [];
            for j = 1:length(CLUSTER_SWEEP_VALS)
                %-
                clust_i = CLUSTER_SWEEP_VALS(j);
                
                %## Calculate dipole positions
                cluster_update = cluster_comp_dipole(TMP_ALLEEG, all_solutions{j}.solutions{1});
                TMP_STUDY.cluster = cluster_update;
                %## (PLOT) Looking for dipoles outside of brain.
                vol = load(MNI_VOL);
                try
                    vol = vol.vol;
                catch
                    vol = vol.mesh;
                end
                rmv_dip = [];
                fig = figure('color','w');
                ft_plot_mesh(vol.bnd(3));
                hold on;
                for c = 2:length(TMP_STUDY.cluster)
                    for d = 1:size(TMP_STUDY.cluster(c).all_diplocs,1)
                        depth = ft_sourcedepth(TMP_STUDY.cluster(c).all_diplocs(d,:), vol);
                        if depth > 20 %give it some leaway 1cm % positive if outside, negative if inside
                            rmv_dip = [rmv_dip;[c,d]];
                            plot3(TMP_STUDY.cluster(c).all_diplocs(d,1),TMP_STUDY.cluster(c).all_diplocs(d,2),TMP_STUDY.cluster(c).all_diplocs(d,3),'*-');
                        end
                    end
                end
                hold off;
                drawnow;
                mkdir([tmp_dir filesep num2str(clust_i)]);
                saveas(fig,[tmp_dir filesep num2str(clust_i) filesep sprintf('ics_out_of_brain.fig')]);
                if ~isempty(rmv_dip)
                    sets_ob = zeros(size(rmv_dip,1),1);
                    comps_ob = zeros(size(rmv_dip,1),1);
                    clusts = unique(rmv_dip(:,1));
                    cnt = 1;
                    for c_i = 1:length(clusts)
                        inds = clusts(c_i)==rmv_dip(:,1);
                        d_i = rmv_dip(inds,2);
                        sets_ob(cnt:cnt+length(d_i)-1) = TMP_STUDY.cluster(clusts(c_i)).sets(d_i);
                        comps_ob(cnt:cnt+length(d_i)-1) = TMP_STUDY.cluster(clusts(c_i)).comps(d_i);
                        TMP_STUDY.cluster(clusts(c_i)).sets(d_i) = [];
                        TMP_STUDY.cluster(clusts(c_i)).comps(d_i) = [];
                        cnt = cnt + length(d_i);
                    end
                    TMP_STUDY.cluster(end+1).sets = sets_ob';
                    TMP_STUDY.cluster(end).comps = comps_ob';
                    TMP_STUDY.cluster(end).name = 'Outlier cluster_outside-brain';
                    TMP_STUDY.cluster(end).parent = TMP_STUDY.cluster(2).parent;
                    TMP_STUDY.cluster(end).algorithm = 'ft_sourcedepth < 0';
                end
                %## REMOVE BASED ON RV
                % [cluster_update] = evaluate_cluster(STUDY,ALLEEG,clustering_solutions,'min_rv');
                cluster_pre_reduce = TMP_STUDY.cluster;
                reduce_method_all = {'min_rv','min_ic','pca_reduce','max_iclabel'};
                
                for p = 1:4 %3 reduce methods
                    reduce_method = reduce_method_all{p}
                    fprintf('Reducing by %s \n',reduce_method);
                    cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method filesep 'in_brain'];
                    if ~exist(cluster_dir,'dir')
                        mkdir(cluster_dir)
                    end
                    TMP_STUDY_TMP = TMP_STUDY;
                    switch reduce_method
                        case 'min_rv'
                            [TMP_STUDY_TMP,~,~] = cluster_rv_reduce(TMP_STUDY_TMP,TMP_ALLEEG);
                        case 'min_ic'
                            [TMP_STUDY_TMP,~,~] = cluster_ica_reduce(TMP_STUDY_TMP,TMP_ALLEEG);
                        case 'pca_reduce'
                            TMP_STUDY_TMP = TMP_STUDY;
                        case 'max_iclabel'
                            [TMP_STUDY_TMP,~,~] = cluster_iclabel_reduce(TMP_STUDY_TMP,TMP_ALLEEG);
                    end
                    cluster_update = TMP_STUDY_TMP.cluster;
                    %- get cluster centroid and residual variance
                    cluster_update = cluster_comp_dipole(TMP_ALLEEG, cluster_update);
                    save([cluster_dir filesep sprintf('cluster_update_spca_%s_%i_%s.mat',study_fName_1,clust_i,reduce_method)],'cluster_update','cluster_pre_reduce');
                    TMP_STUDY_TMP.cluster = cluster_update;

                    %## Look up cluster centroid Brodmann area
                    [~,atlas_names,~] = add_anatomical_labels(TMP_STUDY_TMP,TMP_ALLEEG);
                    for k = 1:length(cluster_update)
                        cluster_update(k).analabel = atlas_names{k,2};
                    end
                    par_save(cluster_update,cluster_dir,sprintf('cl_inf_spca_inbrain_%s_%i_%s.mat',study_fName_1,clust_i,reduce_method));
                    TMP_STUDY_TMP.etc.cluster_vars(j).fpath = [cluster_dir filesep sprintf('cl_inf_spca_inbrain_%s_%i_%s.mat',study_fName_1,clust_i,reduce_method)];
                    TMP_STUDY_TMP.etc.cluster_vars(j).inf = {'type','kmeans','k',clust_i,'params',CLUSTER_PARAMS,'reduce_method',reduce_method};

                %## (PLOT) Looking for dipoles outside of brain.
                %{
                vol = load(MNI_VOL);
                vol = vol.vol;
                vol_bounds = zeros(3,2,length(vol.bnd));
                rmv_dip = [];
                fig = figure('color','w');
        %         ft_plot_mesh(vol.bnd(1));
        %         ft_plot_mesh(vol.bnd(2));
                ft_plot_mesh(vol.bnd(3));
                hold on;
                for c = 2:length(TMP_STUDY.cluster)
                    for d = 1:size(TMP_STUDY.cluster(c).all_diplocs,1)
                        for b = 1:length(vol.bnd)
                            dist = vol.bnd(b).pnt-repmat(TMP_STUDY.cluster(c).all_diplocs(d,:),size(vol.bnd(b).pnt,1),1);
                            sig = sign(dist);
                            pos = sig > 0;
                            neg = sig < 0;
                            depth = ft_sourcedepth(TMP_STUDY.cluster(c).all_diplocs(d,:), vol);
        %                     disp(depth);
                            if depth > 0 && b%any(all(pos,1)) || any(all(neg,1))
                                rmv_dip = [rmv_dip;[c,d,b]];
                                plot3(TMP_STUDY.cluster(c).all_diplocs(d,1),TMP_STUDY.cluster(c).all_diplocs(d,2),TMP_STUDY.cluster(c).all_diplocs(d,3),'*-');
                            end
                        end
                    end
                end
                hold off;
                drawnow;
                saveas(fig,[tmp_dir filesep sprintf('ics_out_of_brain.fig')]);
                %}
                %- save
                    [TMP_STUDY_TMP,TMP_ALLEEG] = parfunc_save_study(TMP_STUDY_TMP,TMP_ALLEEG,...
                                                        TMP_STUDY_TMP.filename,TMP_STUDY_TMP.filepath,...
                                                        'RESAVE_DATASETS','off');
                end
            end
        end
    end
end
