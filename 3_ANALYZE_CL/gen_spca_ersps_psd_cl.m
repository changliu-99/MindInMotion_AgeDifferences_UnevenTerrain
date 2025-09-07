%   Project Title: MiM Age difference paper
%   
%   Perform SPCA on the dataset
%
%   Code Designer: Jacob salminen and Chang Liu
%
%   Version History --> See details at the end of the script.
%   Current Version:  v1.0.20220103.0
%   Previous Version: n/a
%   Summary: 

% sbatch /blue/dferris/jsalminen/GitHub/par_EEGProcessing/src/3_ANALYZE/MIM_OA/run_gen_spca_ersps_timewarp.sh

%{
%## RESTORE MATLAB
% WARNING: restores default pathing to matlab 
restoredefaultpath;
clc;
close all;
clearvars
%}
% Chang Liu - 20240212 run SPCA on the PSD
clear
%% Initialization
% opengl('dsave', 'software') % might be needed to plot dipole plots?
%## TIME
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
colormap_ersp = othercolor('RdYlBu11');
colormap_ersp = colormap_ersp(end:-1:1,:);

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
COND_DESIGNS = {'flat','low','med','high','0p25','0p5','0p75','1p0'};
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
dt = '12012023_OAYA104_icc0p65-0p4_changparams';
%## soft define
DATA_DIR = [source_dir ];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);

study_fName_1 = 'rest_study_new_all_comps_fixed';
% TRIAL_OVERRIDE_FPATH = [STUDIES_DIR filesep 'subject_mgmt' filesep 'trial_event_indices_override.xlsx'];

save_dir = [STUDIES_DIR filesep sprintf('%s',dt) filesep 'cluster'];
load_dir = [STUDIES_DIR filesep sprintf('%s',dt)];
%- create new study directory
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% ===================================================================== %%
%## LOAD STUDIES && ALLEEGS
%- Create STUDY & ALLEEG structs
OVERWRITE_study = true;
if ~exist([load_dir filesep study_fName_1 '.study'],'file')|OVERWRITE_study    
    % Create study - Chang Liu
    [STUDY,ALLEEG] = mim_create_epoch_study_ersp([],'epoch_study_new_all_comps_fixed',load_dir,['GAIT_EPOCHED_ALL' filesep '0p250p50p751p0flatlowmedhigh'],'GAIT_EPOCHED_ALL','SLIDING_EPOCHED_ALL');

%     [STUDY,ALLEEG] = parfunc_save_study(STUDY,ALLEEG,...
%                                             STUDY.filename,STUDY.filepath,...
%                                             'RESAVE_DATASETS','off'); 
    
    [STUDY_REST,ALLEEG_REST] = mim_create_epoch_study_ersp([],study_fName_1,load_dir,['SLIDING_EPOCHED_ALL' filesep 'rest'],'GAIT_EPOCHED_ALL','SLIDING_EPOCHED_ALL');   
%     [STUDY_REST,ALLEEG_REST] = parfunc_save_study(STUDY_REST,ALLEEG_REST,...
%                                             STUDY_REST.filename,STUDY_REST.filepath,...
%                                             'RESAVE_DATASETS','off');  
    if length(ALLEEG) ~= length(ALLEEG_REST)
        fprintf('ALLEEG doesnt match');
        exit()
    end
end

%% ===================================================================== %%

%% CALCULATE GRANDAVERAGE WARPTOs
for subj_i = 1:length(ALLEEG)
    %- assign percondition timewarping
    ALLEEG(subj_i).timewarp.warpto = nanmedian(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
%     ALLEEG(subj_i).timewarp.warpto = nanmean(cat(1,ALLEEG(subj_i).etc.timewarp_by_cond.warpto));
end
allWarpTo = nan(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
% allWarpTo = zeros(length(ALLEEG),size(ALLEEG(1).timewarp.warpto,2));
for subj_i = 1:length(ALLEEG)
    allWarpTo(subj_i,:) = ALLEEG(subj_i).timewarp.warpto; %stack subject specific median event latencies
end
% grandAvgWarpTo = floor(nanmedian(allWarpTo)); % tends to be shorter? (e.g., [0,242,686,915,1357])
% averaged_warpto_events = floor(nanmean(allWarpTo,1)); % tends to be longer? (e.g., [0,262,706,982,1415])
% save([STUDY.filepath filesep 'averaged_warpto_events.mat'],'averaged_warpto_events');
averaged_warpto_events = [0 256 728 999 1461]; % set to be number to ease computation
%% (ERSP PLOT PREP) PREPARE STUDYFILE FOR EXTRACTION (BLACK-HAWK DOWN!)
TIMEWARP_NTIMES = floor(ALLEEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles
ERSP_CROP_TIMES=[averaged_warpto_events(1), averaged_warpto_events(end)+1];
STUDY.etc.averaged_warpto_events = averaged_warpto_events;
fprintf('Using timewarp limits: [%0.4g,%0.4f]\n',averaged_warpto_events(1),averaged_warpto_events(end));
disp(averaged_warpto_events);
%## ersp plot per cluster per condition

%% (PRECOMPUTE MEASURES) COMPUTE ERSPs for Resting conditions
%## ersp plot per cluster per condition
OVERWRITE_SPEC = 0;
DO_BASELINE_CORRECTION = 0;

STUDY_REST = pop_statparams(STUDY_REST,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY_REST = pop_erspparams(STUDY_REST,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY_REST = pop_specparams(STUDY_REST,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');

STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
STUDY = pop_erspparams(STUDY,'subbaseline',ERSP_PARAMS.subbaseline,...
      'ersplim',ERSP_PARAMS.ersplim,'freqrange',ERSP_PARAMS.freqrange,'timerange',ERSP_CROP_TIMES);
SPEC_PARAMS.subtractsubjectmean = 'on';
STUDY = pop_specparams(STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
    'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
    'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');

disp(['Grand average (across all subj) warp to: ',num2str(averaged_warpto_events)]);
DO_TIMEWARP = 1;
DO_BASELINE_CORRECTION = 0;

%%
% parfor (subj_i = 1,SLURM_POOL_SIZE)
for subj_i = 1:length(ALLEEG)
    clear restSPEC gaitSPEC spec 
    EEG = ALLEEG_REST(subj_i);
    icatimf_f = [EEG.filepath filesep sprintf('restSPEC_%s.mat',EEG.subject)];
    
        TMP_STUDY = STUDY_REST;
        EEG = eeg_checkset(EEG,'loaddata');
        if isempty(EEG.icaact)
            fprintf('%s) Recalculating ICA activations\n',EEG.subject);
            EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
            EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
        end

        %- overrride datasetinfo to trick std_precomp to run.
        TMP_STUDY.datasetinfo = STUDY_REST.datasetinfo(subj_i);
        TMP_STUDY.datasetinfo(1).index = 1;
        %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
%             [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
%                     'recompute','on','ersp','on','itc','off',...
%                     'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
%                     'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
%                     'ntimesout',TIMEWARP_NTIMES,...
%                     'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
%             [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
%                     'recompute','on','ersp','on','itc','off',...
%                     'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
%                     'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES}); %ERSP
             filepath = TMP_STUDY.datasetinfo(1).filepath;
             trialinfo = std_combtrialinfo(TMP_STUDY.datasetinfo, 1);
             filebase = fullfile(filepath, EEG.subject);

%               [~, times,logfreqs,parameters,logersp] = std_ersp_cl(EEG,'components',1:size(EEG.icawinv,2),'savetrials','off',...
%                       'recompute','on','fileout', filebase, 'trialinfo', trialinfo,...
%                       'parallel','on','cycles',ERSP_PARAMS.cycles,...
%                       'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
%                       'baseline',nan());%logersp is frequency x time x component
            [spec, freqs] = std_spec(EEG, 'components',1:size(EEG.icawinv,2),...
                'savetrials','on','recompute','on','fileout', filebase, 'trialinfo', trialinfo,...
                'specmode','psd','logtrials','off','subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
                'freqrange',SPEC_PARAMS.plot_freqrange,'savefile','off');
        end
        fprintf('Done calculating timewarped ERSPs for resting condition %s \n',EEG.subject);
        
        restSPEC.freqs = freqs;        
        restSPEC.spec(1,:,:) = mean(spec,3);% time x comps x freq
        restSPEC.logspec(1,:,:) = 10*log10(mean(spec,3));
        
        save([EEG.filepath filesep ['restSPEC_',EEG.subject,'.mat']],'restSPEC');
        
        figure();plot(freqs,squeeze(restSPEC.logspec(1,:,:)));
%     [allersp_rest, alltimes_rest, allfreqs_rest] = compute_ersp(STUDY_REST,EEG);%allersp is freq x time x channel
%     allersp_rest_avg = squeeze(mean(allersp_rest{1},2));


%% (PRECOMPUTE MEASURES) COMPUTE ERSPs
% parfor (subj_i = 1,SLURM_POOL_SIZE)
    
    EEG = ALLEEG(subj_i);
    icatimf_f = [EEG.filepath filesep sprintf('gaitSPEC_%s.mat',EEG.subject)];
    
    TMP_STUDY = STUDY;
    EEG = eeg_checkset(EEG,'loaddata');
    if isempty(EEG.icaact)
        fprintf('%s) Recalculating ICA activations\n',EEG.subject);
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end

    %- overrride datasetinfo to trick std_precomp to run.
    TMP_STUDY.datasetinfo = STUDY.datasetinfo(subj_i);
    TMP_STUDY.datasetinfo(1).index = 1;
    %- determine timewarping parameters
     if DO_TIMEWARP
        timewarp_param = EEG.timewarp.latencies;
        timewarpms_param = averaged_warpto_events;
     else
        timewarp_param = [];
        timewarpms_param = [];
     end
     if ~exist(icatimf_f,'file') | OVERWRITE_SPEC% || any(strcmp(EEG.subject,FINISHED_ADULTS))
         %-
        if DO_BASELINE_CORRECTION
            % Baseline correction
%             [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
%                     'recompute','on','ersp','on','itc','off',...
%                     'erspparams',{'parallel','off','cycles',ERSP_PARAMS.cycles,...
%                     'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),...
%                     'ntimesout',TIMEWARP_NTIMES,'timewarp',timewarp_param,...
%                     'timewarpms',timewarpms_param,'baseline',[averaged_warpto_events(1),averaged_warpto_events(end)],...
%                     'trialbase','off','basenorm','on'}); %ERSP
        else
            % No baseline correction
%             [~, ~] = std_precomp(TMP_STUDY,EEG,'components','savetrials','on',...
%                     'recompute','on','ersp','on','itc','off',...
%                     'erspparams',{'parallel','on','cycles',ERSP_PARAMS.cycles,...
%                     'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
%                     'baseline',nan(),'timewarp',timewarp_param,...
%                     'timewarpms',timewarpms_param}); %ERSP
            filepath = TMP_STUDY.datasetinfo(1).filepath;
            trialinfo = std_combtrialinfo(TMP_STUDY.datasetinfo, 1);
            filebase = fullfile(filepath, EEG.subject);

%             EEG = pop_selectevent(EEG, 'cond','high','deleteevents','off','deleteepochs','on','invertepochs','off'); 1:size(EEG.icawinv,2)
%             [~, times,logfreqs,parameters,logersp] = std_ersp_cl(EEG,'components',1:size(EEG.icawinv,2),'savetrials','off',...
%                       'recompute','on','fileout', filebase, 'trialinfo', trialinfo,...
%                       'parallel','on','cycles',ERSP_PARAMS.cycles,...
%                       'nfreqs',length((ERSP_PARAMS.freqrange(1):ERSP_PARAMS.freqrange(2))),'ntimesout',TIMEWARP_NTIMES,...
%                       'baseline',nan(),'timewarp',timewarp_param,...
%                       'timewarpms',timewarpms_param);%logersp is frequency x time x component
            
            [spec, freqs] = std_spec(EEG, 'components',1:size(EEG.icawinv,2),...
                'savetrials','on','recompute','on','fileout', filebase, 'trialinfo', trialinfo,...
                'specmode','psd','logtrials','off','subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
                'freqrange',SPEC_PARAMS.plot_freqrange,'savefile','off');
            
            gaitSPEC.freqs = freqs;        
            gaitSPEC.spec(1,:,:) = mean(spec,3);% time x comps x freq
            gaitSPEC.logspec(1,:,:) = 10*log10(mean(spec,3));
        end
        fprintf('Done calculating timewarped ERSPs for %s',EEG.subject);
     
        %% subtract baseline
        gaitSPEC_subRest = bsxfun(@minus, gaitSPEC.logspec, restSPEC.logspec(:,:,:)); %spec is 1 x comps x freq
        %required by spca times x component x freq
%               
        % version - do not subtract rest
        %% get the denoise coefficient
        [SPEC_corr, ~, PSC1, ~,V] = specPCAdenoising_CL(gaitSPEC_subRest);

        gaitSPEC.gaitSPEC_subRest = gaitSPEC_subRest;
        gaitSPEC.SPEC_corr = SPEC_corr;
        gaitSPEC.GPM_corr = [];
        gaitSPEC.PSC1 = PSC1;
        gaitSPEC.V = V;
        gaitSPEC.GPM = [];
        
        save([EEG.filepath filesep ['gaitSPEC_',EEG.subject,'.mat']],'gaitSPEC');
    else
        fprintf('Timewarped SPEC already calculated for %s \n',EEG.subject);
        load([EEG.filepath filesep ['gaitSPEC_',EEG.subject,'.mat']]);
    end
%     [allersp_gait, alltimes_gait, allfreqs_gait] = compute_ersp(STUDY,EEG);%allersp is freq x time x channel
%     save([EEG.filepath filesep 'allersp_gait.mat'],'allersp_gait','alltimes_gait');
    
    %% Apply this coefficient V to each condition
    
    for ci = 1:length(COND_DESIGNS)
        clear gaitSPEC_cond
        try
            % clean each condition with spca
            fprintf('=== Cleaning %s for %s ====\n',COND_DESIGNS{ci},EEG.subject);

            tempEEG = pop_selectevent(EEG, 'cond',COND_DESIGNS{ci},'deleteevents','off','deleteepochs','on','invertepochs','off');
            timewarp = make_timewarp(tempEEG,{'RHS','LTO','LHS','RTO','RHS'},'baselineLatency',0, ...
            'maxSTDForAbsolute',Inf,...
            'maxSTDForRelative',Inf);
        
            tempEEG.timewarp = timewarp;
            % important - forget to timewarp
            timewarp_param = tempEEG.timewarp.latencies;
            timewarpms_param = averaged_warpto_events;
            
            filepath = TMP_STUDY.datasetinfo(1).filepath;
            trialinfo = std_combtrialinfo(TMP_STUDY.datasetinfo, 1);
            filebase = fullfile(filepath, tempEEG.subject);
            [spec, freqs] = std_spec(tempEEG, 'components',1:size(EEG.icawinv,2),...
                'savetrials','on','recompute','on','fileout', filebase, 'trialinfo', trialinfo,...
                'specmode','psd','logtrials','off','subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
                'freqrange',SPEC_PARAMS.plot_freqrange,'savefile','off');          

            gaitSPEC_cond.freqs = freqs;        
            gaitSPEC_cond.spec(1,:,:) = mean(spec,3);% time x comps x freq
            gaitSPEC_cond.logspec(1,:,:) = 10*log10(mean(spec,3));

            gaitSPEC_subRest = bsxfun(@minus, gaitSPEC_cond.logspec, restSPEC.logspec(:,:,:));
%             gaitSPEC_subRest = permute(gaitSPEC_subRest,[2,3,1]);%required by spca times x component x freq
%             GPM = bsxfun(@minus,gaitSPEC_subRest,mean(gaitSPEC_subRest));        

            [SPEC_corr, GPM_corr, PSC1] = specPCAdenoising_CL(gaitSPEC_subRest, gaitSPEC.V);% V matrix:  freqs x freqs 

            gaitSPEC_cond.ID         = EEG.subject;
            gaitSPEC_cond.cond       = COND_DESIGNS{ci};
            gaitSPEC_cond.times      = 1;
            gaitSPEC_cond.freqs      = freqs;  
            
            gaitSPEC_cond.gaitSPEC_subRest    = gaitSPEC_subRest;
            gaitSPEC_cond.SPEC_corr  = SPEC_corr;
            gaitSPEC_cond.GPM_corr   = [];
            gaitSPEC_cond.PSC1       = PSC1;
            gaitSPEC_cond.GPM        = [];

            save([EEG.filepath filesep ['gaitSPEC_',EEG.subject,'_',COND_DESIGNS{ci},'.mat']],'gaitSPEC_cond');

            %% plot figures
            % ERSP uncorrected
            cfg.f_axis = gaitSPEC_cond.freqs;
            cfg.t_axis = gaitSPEC_cond.times;
            clim = 9;
            f1 = figure();subplot(1,2,1);
            data = squeeze(mean(gaitSPEC_subRest,2));
            plot(cfg.f_axis,data'); hold on % plot ERSP
            xlim([3 100]);
            
            subplot(1,2,2);
            data = squeeze(mean(SPEC_corr,2));
            plot(cfg.f_axis,data'); hold on % plot ERSP  
            xlim([3 100]);
            saveas(f1,[EEG.filepath filesep ['SPEC_',EEG.subject,'_',COND_DESIGNS{ci},'_figure.jpg']])

 
        catch e
            fprintf('Cannot process condition %s \n',COND_DESIGNS{ci});
            fprintf('%s\n',getReport(e))
        end
        close all;
    end
end

                     
