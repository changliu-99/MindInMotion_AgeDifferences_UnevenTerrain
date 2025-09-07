% Script to gather PSD, ERSP outcome after SPCA and clustering 
% Run fooof analysis on PSDs after cleaned with sPCA

% Chang Liu - 20240311

clear all;close all

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
% cd(run_dir)
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
CLUSTER_SWEEP_VALS = [9,10,11,12,13];%[12,14,16,19]; %[10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
% DO_K_DISTPRUNE = false;
DO_K_ICPRUNE = true;
% DO_K_SWEEPING = false;
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
dt = '12012023_OAYA104_icc0p65-0p4_changparams';
%## soft define
DATA_DIR = [source_dir ];
STUDIES_DIR = [DATA_DIR filesep DATA_SET filesep '_studies'];
% study_fName_1 = sprintf('%s_EPOCH_study',[TRIAL_TYPES{:}]);

study_fName_1 = 'rest_study_new_spca_all_fixed_reduce';
study_fName_1_all = 'gait_study_new_spca_all_fixed_reduce';
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
% if ~ispc
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
% else
%     [STUDY,ALLEEG] = pop_loadstudy('filename',[study_fName_1 '.study'],'filepath',load_dir);
% end
% 
% CLUSTER_PARAMS.filename = STUDY.filename;
% CLUSTER_PARAMS.filepath = STUDY.filepath;

%% Use the rest study  
study_fName = sprintf('temp_study_rejics_spca_%s_%i',study_fName_1, MIN_ICS_SUBJ);
tmp_dir = [save_dir filesep sprintf('icrej_spca_%s_%i',study_fName_1,MIN_ICS_SUBJ)];
if ~exist([tmp_dir filesep study_fName '.study'],'file')
    error('ERROR. study file does not exist');
else
    if ~ispc
        [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '_UNIX.study'],'filepath',tmp_dir);
    else
        TMP_STUDY = importdata(fullfile(tmp_dir,[study_fName '.study']));
        TMP_STUDY.filepath = convertPath2Drive(TMP_STUDY.filepath);
        for i = 1:length(TMP_STUDY.datasetinfo)
            TMP_STUDY.datasetinfo(i).filepath = convertPath2Drive(TMP_STUDY.datasetinfo(i).filepath);
        end
        %%%!!!!^^%%%% Superimportant 
        % somehow the TMP_STUDY.subject does not match with
        % TMP_STUDY.datasetinfo.subject - 2025.05.07 Chang Liu, make manual
        % correction here
        TMP_STUDY.subject = {TMP_STUDY.datasetinfo.subject};
        
        
        STUDY = TMP_STUDY;
        if ispc
            save([tmp_dir filesep sprintf('temp_study_rejics_spca_%s_%i.study',study_fName_1, MIN_ICS_SUBJ)],'STUDY');
        else
            save([tmp_dir filesep sprintf('temp_study_rejics_spca_%s_%i_UNIX.study',study_fName_1, MIN_ICS_SUBJ)],'STUDY');
        end
        [TMP_STUDY,TMP_ALLEEG] = pop_loadstudy('filename',[study_fName '.study'],'filepath',tmp_dir);
        
        TMP_STUDY.subject = {TMP_STUDY.datasetinfo.subject};
    end
end
STUDY_GROUP_DESI = {{'subjselect',{},...
    'variable1','cond','values1',{'flat','low','med','high'},...
    'variable2','group','values2',unique({TMP_STUDY.datasetinfo.group})},...
    {'subjselect',{},...
    'variable2','cond','values1',{'0p25','0p5','0p75','1p0'},...
    'variable2','group','values2',unique({TMP_STUDY.datasetinfo.group})}};
%## grab subjects for study designs
tmp_group_orig = cell(length(TMP_ALLEEG),1);
tmp_group_unif = cell(length(TMP_ALLEEG),1);
for subj_i = 1:length(TMP_ALLEEG)
    tmp_group_orig{subj_i} = TMP_ALLEEG(subj_i).group;
    tmp_group_unif{subj_i} = 'Older Adults';
end
%## ersp plot per cluster per condition
TMP_STUDY = pop_statparams(TMP_STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);
TMP_STUDY = pop_specparams(TMP_STUDY,'subtractsubjectmean',SPEC_PARAMS.subtractsubjectmean,...
        'freqrange',SPEC_PARAMS.plot_freqrange,'plotmode','condensed',...
        'plotconditions','together','ylim',SPEC_PARAMS.plot_ylim,'plotgroups','together');
    
%%%%%% Transfer information to walking trial study
if ispc
    STUDY_gait = importdata(fullfile(load_dir,['epoch_study_spca_all','.study']));
    STUDY_gait.filepath = convertPath2Drive(STUDY_gait.filepath);
    
    STUDY_gait.subject = {STUDY_gait.datasetinfo.subject};
else
    STUDY_gait = importdata(fullfile(load_dir,['epoch_study_spca_all_UNIX','.study']));
    STUDY_gait.subject = {STUDY_gait.datasetinfo.subject};
end
TMP_STUDY_gait = TMP_STUDY;
% replace subject information 
for n = 1:length(TMP_STUDY_gait.datasetinfo)
    n_ind = find(strcmp({STUDY_gait.datasetinfo(:).subject},TMP_STUDY.datasetinfo(n).subject));
    TMP_STUDY_gait.datasetinfo(n).filepath = STUDY_gait.datasetinfo(n_ind).filepath;
    TMP_STUDY_gait.datasetinfo(n).filename = STUDY_gait.datasetinfo(n_ind).filename;
    TMP_STUDY_gait.datasetinfo(n).condition = 'gait';    
    TMP_STUDY_gait.datasetinfo(n).trialinfo = STUDY_gait.datasetinfo(n_ind).trialinfo;
end
TMP_STUDY_gait.name = sprintf('temp_study_rejics_spca_%s_%i',study_fName_1_all, MIN_ICS_SUBJ);
STUDY = TMP_STUDY_gait;
if ispc
    save([tmp_dir filesep sprintf('temp_study_rejics_spca_%s_%i.study',study_fName_1_all, MIN_ICS_SUBJ)],'STUDY');
else
    save([tmp_dir filesep sprintf('temp_study_rejics_spca_%s_%i_UNIX.study',study_fName_1_all, MIN_ICS_SUBJ)],'STUDY');
end
%%%%%
[TMP_STUDY_gait,TMP_ALLEEG_gait] = pop_loadstudy('filename',[sprintf('temp_study_rejics_spca_%s_%i.study',study_fName_1_all, MIN_ICS_SUBJ)],'filepath',tmp_dir);
TMP_STUDY_gait.subject = {TMP_STUDY_gait.datasetinfo.subject};


TMP_STUDY_gait.design = [];
TMP_STUDY_gait = std_makedesign(TMP_STUDY_gait, TMP_ALLEEG_gait, 1, 'subjselect',{},'variable1','cond','values1', {'flat','low','med','high'});
TMP_STUDY_gait.design(1).variable(2) = [];
%--- Add design information
TMP_STUDY.design = [];
TMP_STUDY = std_makedesign(TMP_STUDY, TMP_ALLEEG, 1,'variable1','cond','values1', {'rest'});

% TMP_STUDY_gait = std_makedesign(TMP_STUDY_gait, TMP_ALLEEG_gait, 2, 'variable1','cond','values1', {'0p25','0p5','0p75','1p0'});

% TMP_STUDY_gait = std_makedesign(STUDY, ALLEEG, 2, 'variable1','cond','values1', {'0p25','0p5','0p75','1p0'},'variable2','group','values1',unique({TMP_STUDY_gait.datasetinfo.group}));

%% Compute scalp on
[~, ~] = std_precomp(TMP_STUDY, TMP_ALLEEG,...
                        'components',...
                        'spec','off',...
                        'scalp','on',...
                        'recompute','off');

%% Check dataset to see if it has all conditions
condition_store = cell(length(TMP_ALLEEG_gait),8);
for i = 1:length(TMP_ALLEEG_gait)
    condition_store(i,1:length(unique({TMP_ALLEEG_gait(i).event.cond}))) = unique({TMP_ALLEEG_gait(i).event.cond});
end

%% Generate the cluster topography
reduce_method = 'max_iclabel';
if isunix;addpath '/blue/dferris/liu.chang1/scripts/MiM_HY/_common/std_plot_funs/';end

for j = 3%:length(CLUSTER_SWEEP_VALS)% cluster = 12
    %-
    fprintf('Generating the cluster = %s \n',num2str(CLUSTER_SWEEP_VALS(j)));
    clust_i = CLUSTER_SWEEP_VALS(j);
    cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];
    cluster_update = load([cluster_dir filesep 'in_brain' filesep sprintf('cl_inf_spca_inbrain_%s_%i_%s.mat',study_fName_1,clust_i,reduce_method)]); %par_load(cluster_dir,sprintf('cluster_inf_%i.mat',clust_i));
    
    %## Create Plots
    TMP_STUDY.cluster = cluster_update.SAVEVAR;
    TMP_STUDY_gait.cluster = cluster_update.SAVEVAR;
    
    ics_subjs = cellfun(@(x) size(x,2),{TMP_STUDY.datasetinfo.comps});
    [h,p,~,stats] = ttest2(ics_subjs(1:31),ics_subjs(32:end));
    
     %- get inds
    [~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps_CL(TMP_STUDY);
%     mkdir(fullfile(cluster_dir,'ERSP_Plots'));
%     save(fullfile(cluster_dir,'ERSP_Plots','cluster_info.mat'),'valid_clusters','main_cl_inds');
    
    %- clusters to plot
    CLUSTER_PICKS = main_cl_inds(2:end); 
    %## PLOT CLUSTER BASE INFORMATION
    %- CLUSTER DIPOLES, TOPOS
    fprintf('There are a total of %i young adults and %i high functioning older adults and %i frail older adults \n',...
        sum(contains({TMP_STUDY.datasetinfo(:).group},'H1')),...
        sum(contains({TMP_STUDY.datasetinfo(:).group},'H2')),...
        sum(contains({TMP_STUDY.datasetinfo(:).group},'H3')));

    for file = 1:length(TMP_ALLEEG)
        if isunix
            TMP_ALLEEG(file).dipfit.mrifile = '/blue/dferris/liu.chang1/scripts/MiM_CRUNCH/_resources/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii';
        elseif ispc
            TMP_ALLEEG(file).dipfit.mrifile = 'M:\liu.chang1\scripts\MiM_CRUNCH\_resources\mni_icbm152_nlin_sym_09a\mni_icbm152_t1_tal_nlin_sym_09a.nii';
        end
    end
    mim_gen_cluster_figs_CL(TMP_STUDY,TMP_ALLEEG,cluster_dir,...
        'CLUSTERS_TO_PLOT', valid_clusters,...
        'DO_SINGLE_CLUSTER_PLOTS',1,...
        'DO_SINGLE_CLUSTER_PLOTS_Group',1,...
        'DO_GROUP_DIPOLE',1);
%     
%     % get the aggregate anatomical label for one cluster
%     [atlas_names] = aggregate_anatomical_labels(cluster_update.SAVEVAR,16);
%     [uni,~,idx] = unique(cellstr(char(atlas_names{:,2})));
%     hist(idx,unique(idx))
%     xticklabels(uni); 
%     xtickangle(45)
    %- close all figures
    close all
end

% CLUSTER_SELECT = input('What cluster do you choose for analysis?  ');
CLUSTER_SELECT = CLUSTER_SWEEP_VALS(j);

fprintf('Processing all for %s \n',num2str(CLUSTER_SELECT));
%% ====== Gather ERSP for each cluster
%% Gather ERSP for each cluster
gatherERSP = 1;
if gatherERSP
    group_list = {TMP_STUDY.datasetinfo(:).group};
    for i = valid_clusters
        tic
        fprintf('\n grabbing ersp for %s \n',num2str(i));
        sub_cl = TMP_STUDY.cluster(i).sets;
        sub_comp = TMP_STUDY.cluster(i).comps;
        j = 1;
        ERSP_GPM_corr = cell(4,1);
        for n = 1:length(sub_cl)
            for cond = 1:length(terrain_trials)-1
                sub_ID = TMP_STUDY.datasetinfo(sub_cl(n)).subject;%use the datasetinfo instead
                fprintf('\n loading from %s',sub_ID);
                load_ersp_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitERSP_',sub_ID,'_',terrain_trials{cond}]); 
                load(load_ersp_spca);
                [ERSP_GPM_corr{cond,1}(:,j,:)] = gaitERSP_cond.GPM_corr(:,sub_comp(n),:);%ERSP_spca = time x  comp x freqs
                fprintf('. %s',num2str(sub_comp(n)));
                [ERSP_corr{cond,1}(:,j,:)] = gaitERSP_cond.ERSP_corr(:,sub_comp(n),:);
                [ERSP_GPM{cond,1}(:,j,:)] = gaitERSP_cond.GPM(:,sub_comp(n),:);
            end
            j = j + 1;
        end
        TMP_STUDY.cluster(i).group = group_list(TMP_STUDY.cluster(i).sets);

        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H1000s')) = {1};
        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H2000s')) = {2};
        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H3000s')) = {2};
        TMP_STUDY.cluster(i).group = cellfun(@(x) [[]; x], TMP_STUDY.cluster(i).group);

        group_store = TMP_STUDY.cluster(i).group;

        SPCA_results.ERSP_GPM_corr = ERSP_GPM_corr;
        SPCA_results.ERSP_corr = ERSP_corr;
        SPCA_results.ERSP_GPM = ERSP_GPM;
        SPCA_results.cluster = i;
        SPCA_results.group_store = group_store;

        clust_i = CLUSTER_SELECT;
        cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];

        if ~exist(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)]))
            mkdir(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)]))
        end
        save(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)],'SPCA_results.mat'),'SPCA_results');
        clear SPCA_results ERSP_GPM ERSP_corr ERSP_GPM_corr
        toc
    end
% ------- plot GPM corr and non-corr GPM
load(fullfile(load_dir,'H1004','GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh','H1004.icatimef'),'-mat','parameters');
warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
group = 2;
close all;
for k = 4%[3:14]
    f1 = figure('color','w','position',[100 100 1200 600]);
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.mat'));
    for cond = 1:length(terrain_trials)-1
        subplot(2,4,cond)
        tftopo(squeeze(mean(SPCA_results.ERSP_GPM_corr{cond,1}(:,SPCA_results.group_store==group,:),2)),gaitERSP_cond.times,gaitERSP_cond.logfreqs,...
            'title',terrain_trials{cond},'limits',[0 1461 nan nan -0.5 0.5]);
        xlim([0 1461]);

        subplot(2,4,cond+4)
        tftopo(squeeze(mean(SPCA_results.ERSP_GPM{cond,1}(:,SPCA_results.group_store==group,:),2)),gaitERSP_cond.times,gaitERSP_cond.logfreqs,...
            'title',[terrain_trials{cond},'_raw'],'limits',[0 1461 nan nan -0.5 0.5]);
        xlim([0 1461]);
    end
%     savefig(f1, fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.fig'));
    
    for cond = 1:length(terrain_trials)-1
        subplot(2,4,cond)
        tftopo(squeeze(mean(SPCA_results.ERSP_GPM_corr{cond,1}(:,SPCA_results.group_store==group,:),2)),gaitERSP_cond.times,gaitERSP_cond.logfreqs,...
            'title',terrain_trials{cond},'limits',[warpingvalues(1) warpingvalues(end) nan nan nan nan],...
                'vert',warpingvalues(1:5),'logfreq','native');
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
            title(strcat({'Cluster '},num2str(k)));
    %         ylim(log([4 50]))
            cbar('vert');
%         xlim([0 1461]);ylim([3 100]););
               
        subplot(2,4,cond+4)
        tftopo(squeeze(mean(SPCA_results.ERSP_GPM{cond,1}(:,SPCA_results.group_store== group,:),2)),gaitERSP_cond.times,gaitERSP_cond.logfreqs,...
            'title',[terrain_trials{cond},'_raw'],'limits',[warpingvalues(1) warpingvalues(end) nan nan nan nan],...
                'vert',warpingvalues(1:5),'logfreq','native');
            set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency (Hz)'),'fontsize',16,'fontweight','bold');
            xlabel('Time (ms)','Fontsize',16,'fontweight','bold');
            title(strcat({'Cluster '},num2str(k)));
    %         ylim(log([4 50]))
            cbar('vert');
%         xlim([0 1461]);ylim([3 100]););
%         xlim([0 1461]);ylim([3 100]);
    end
%     savefig(f1, fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results_wind.fig'));
end
    %% Compute P2P ERSP
    config.warpingvalues = warpingvalues;
    config.colormap_ersp = colormap_ersp;
    config.title_keyword = {'flat','low','med','high'};
    config.STUDY = STUDY;
    config.runPairwise = 0;
    config.YlimMax = 50;
    config.CLUSTER_SWEEP_VALS = CLUSTER_SELECT;
    config.color = color;

    p2p_ERSP_table_sum = table;
    fprintf('Clusters = %s \n',num2str(config.CLUSTER_SWEEP_VALS));
    for k = valid_clusters

        config.k = k;
        %% LOAD SPCA result
        load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.mat'));
        sub_ID = TMP_STUDY.datasetinfo(1).subject;%use the datasetinfo instead
        load_ersp_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitERSP_',sub_ID,'_','flat']);
        load(load_ersp_spca);
        SPCA_results.times = gaitERSP_cond.times;
        SPCA_results.logfreqs = gaitERSP_cond.logfreqs;

        % First step, run stats for the GPM style data
        allersp = SPCA_results.ERSP_GPM_corr;
        allersp_ERSP = SPCA_results.ERSP_corr;
        for j = 1:4
            allersp_reorder{j,1} = permute(allersp{j},[3 1 2]);%transform to freq x time x subject
            allersp_ERSP_reorder{j,1} = permute(allersp_ERSP{j},[3 1 2]);
        end

        % - Compute the peak to peak ERSP
%         config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'p2p-ERSP.jpg');
%         config.save_file = true;
%         p2p_ERSP_table = compute_PtP_ERSP(allersp_reorder,SPCA_results.times,SPCA_results.logfreqs,config);
%         subj = {TMP_STUDY.datasetinfo(:).subject};
%         p2p_ERSP_table_all = [table(repmat(TMP_STUDY.cluster(k).sets',16,1),'VariableNames',{'subID'})... 
%             table(repmat(subj(TMP_STUDY.cluster(k).sets)',16,1),'VariableNames',{'subjectName'}) ...
%             table(repmat(k,height(p2p_ERSP_table),1),'VariableNames',{'cluster'}) ...
%             table(repmat(SPCA_results.group_store',16,1),'VariableNames',{'group'}) p2p_ERSP_table];
% 
%         p2p_ERSP_table_sum = [p2p_ERSP_table_sum;p2p_ERSP_table_all];

    end
%     p2p_ERSP_table_sum_unstack = unstack(p2p_ERSP_table_sum,{'p2p_ersp'},'band','VariableNamingRule','preserve');
% 
%     % -- SAVE simplified table that only contained relevant variables - update 20240425
%     spec_data_dir = [cluster_dir filesep 'in_brain' filesep 'spec_data'];
%     save_dir = [spec_data_dir filesep 'psd_calcs_spca'];
%     mkdir(save_dir)
%     
%     save_filename.xlsx_name = 'p2p_ERSP_table_sum';
%     save_filename.psd_table_name = 'p2p_ERSP_table_sum';
%     save_filename.save_dir = save_dir;
% 
%     
%     writetable(p2p_ERSP_table_sum,[save_filename.save_dir filesep [save_filename.xlsx_name,'.xlsx']]);
%     save([save_filename.save_dir filesep [save_filename.psd_table_name,'.mat']],'p2p_ERSP_table_sum');
% 
%     writetable(p2p_ERSP_table_sum_unstack,[save_filename.save_dir filesep [save_filename.xlsx_name,'_unstack.xlsx']]);
%     save([save_filename.save_dir filesep [save_filename.psd_table_name,'_unstack.mat']],'p2p_ERSP_table_sum');

    
end
%% ========= Gather SPEC information
%% Sanity check for SPCA 
% For one example participant, 
plot_supplementary_material = 1;
if plot_supplementary_material
    sub_ID = 'H1007'; % H1004
    j = 1;
    load_spec = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitSPEC_',sub_ID,'.mat']); 
    load(load_spec);
    load_rest = fullfile(load_dir, sub_ID,'SLIDING_EPOCHED_ALL','rest',['restSPEC_',sub_ID,'.mat']);
    load(load_rest);
    SPEC_GPM_corr = cell(4,1);
    % load the IC label result
    load(fullfile(load_dir,sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh','reject_ic',[sub_ID,'_ICRej.mat']));
    muscle_ic = reject_struct.IC_all_muscle >= 2;
    brain_ic = reject_struct.IC_all_brain >= 8;
    % plot the average cleaning result of muscle component
    figure('color','white');
    subplot(1,2,1);hold on;
    title('Muscle')
    x = 0:250;
    y = double(mean(squeeze(gaitSPEC.logspec(1,muscle_ic,:))));
    E = double(std(squeeze(gaitSPEC.logspec(1,muscle_ic,:))));
    JackKnife_sung(x,y',E',color.Blue,color.lightBlue);
    y2 = double(mean(squeeze(gaitSPEC.gaitSPEC_subRest(1,muscle_ic,:))));
    E2 = double(std(squeeze(gaitSPEC.gaitSPEC_subRest(1,muscle_ic,:))));
    JackKnife_sung(x,y2',E2',color.Red,color.lightRed);
    y2 = double(mean(squeeze(gaitSPEC.SPEC_corr(1,muscle_ic,:))));
    E2 = double(std(squeeze(gaitSPEC.SPEC_corr(1,muscle_ic,:))));
    JackKnife_sung(x,y2',E2',color.Green,color.lightGreen);
    y3 = double(mean(squeeze(gaitSPEC.SPEC_corr(1,muscle_ic,:))+squeeze(restSPEC.logspec(1,muscle_ic,:))));
    E3 = double(std(squeeze(gaitSPEC.SPEC_corr(1,muscle_ic,:))+squeeze(restSPEC.logspec(1,muscle_ic,:))));
    JackKnife_sung(x,y3',E3',color.Orange,color.Yellow);
    legend('','Raw','','remove rest','','corrected','','corrected+rest')
    xlim([3 250]);
    xlabel('Frequency(Hz)');ylabel('log Power (dB)');
    set(gca,'fontsize',12);

    subplot(1,2,2);hold on;
    title('Brain')
    x = 0:250;
    y = double(mean(squeeze(gaitSPEC.logspec(1,brain_ic,:))));
    E = double(std(squeeze(gaitSPEC.logspec(1,brain_ic,:))));
    JackKnife_sung(x,y',E',color.Blue,color.lightBlue);
    y2 = double(mean(squeeze(gaitSPEC.gaitSPEC_subRest(1,brain_ic,:))));
    E2 = double(std(squeeze(gaitSPEC.gaitSPEC_subRest(1,brain_ic,:))));
    JackKnife_sung(x,y2',E2',color.Red,color.lightRed);
    y2 = double(mean(squeeze(gaitSPEC.SPEC_corr(1,brain_ic,:))));
    E2 = double(std(squeeze(gaitSPEC.SPEC_corr(1,brain_ic,:))));
    JackKnife_sung(x,y2',E2',color.Green,color.lightGreen);
    y3 = double(mean(squeeze(gaitSPEC.SPEC_corr(1,brain_ic,:))+squeeze(restSPEC.logspec(1,brain_ic,:))));
    E3 = double(std(squeeze(gaitSPEC.SPEC_corr(1,brain_ic,:))+squeeze(restSPEC.logspec(1,brain_ic,:))));
    JackKnife_sung(x,y3',E3',color.Orange,color.Yellow);
    legend('','Raw','','remove rest','','corrected','','corrected+rest')
    xlabel('Frequency(Hz)');ylabel('log Power (dB)');
    xlim([3 250]);
    set(gca,'fontsize',12);

    % plot all brain ics
    figure();
    plot(x,squeeze(gaitSPEC.gaitSPEC_subRest(1,:,:)))

    figure();
    bar(x,gaitSPEC.V(1,:))

    figure();
    plot(x,squeeze(gaitSPEC.PSC1(1,:,:)))

    load_spec_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitSPEC_',sub_ID,'_','med']); 
    load(load_spec_spca);
    [SPEC_corr, GPM_corr, PSC1,allPSC] = specPCAdenoising_CL(gaitSPEC_cond.gaitSPEC_subRest, gaitSPEC.V);

    figure();
    subplot(1,2,1)
    plot(x,squeeze(allPSC{1}(1,brain_ic,:)));
    subplot(1,2,2)
    plot(x,squeeze(allPSC{2}));

    figure()
    imagesc(x,x,squeeze(gaitSPEC.V(:,1,:)))
end
%% Gather SPEC
gatherSPEC = 1;
if gatherSPEC
    all_cond_trials = {'flat','low','med','high','0p25','0p5','0p75','1p0'};
    group_list = {TMP_STUDY.datasetinfo(:).group};
    for i = valid_clusters
        clear SPCA_results_spec SPEC_corr SPEC_corr_addRest SPEC_raw PSC1 SPEC_subRest SPEC_rest
        fprintf('\n grabbing ersp for %s \n',num2str(i));
        sub_cl = TMP_STUDY.cluster(i).sets;
        sub_comp = TMP_STUDY.cluster(i).comps;
        j = 1;
        for n = 1:length(sub_cl)
            sub_ID = TMP_STUDY.datasetinfo(sub_cl(n)).subject;
            fprintf('\n subject no. %s \n component no.%s\n',sub_ID,num2str(sub_comp(n)));
            load_spec = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitSPEC_',sub_ID,'.mat']); 
            load(load_spec);
            load_rest = fullfile(load_dir, sub_ID,'SLIDING_EPOCHED_ALL','rest',['restSPEC_',sub_ID,'.mat']);
            load(load_rest);

            [SPEC_rest{1,1}(:,j,:)] = restSPEC.logspec(:,sub_comp(n),:);
            SPEC_GPM_corr = cell(4,1);  

            for cond = 1:length(all_cond_trials)
                fprintf(all_cond_trials{cond});
                try
                    load_spec_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitSPEC_',sub_ID,'_',all_cond_trials{cond}]); 
                    load(load_spec_spca);

                    %ERSP_spca = time x comp x freqs
                    fprintf('.');
                    [SPEC_corr{cond,1}(:,j,:)] = gaitSPEC_cond.SPEC_corr(:,sub_comp(n),:);

        %             figure();plot(squeeze(gaitSPEC_cond.SPEC_corr(:,sub_comp(n),:)),'r');hold on;plot(squeeze(gaitSPEC_cond.gaitSPEC_subRest(:,sub_comp(n),:)),'b');legend('corrected','raw');
                    [SPEC_corr_addRest{cond,1}(:,j,:)] = gaitSPEC_cond.SPEC_corr(:,sub_comp(n),:) + restSPEC.logspec(:,sub_comp(n),:);

                    [SPEC_raw{cond,1}(:,j,:)] = gaitSPEC_cond.logspec(:,sub_comp(n),:);

                    [PSC1{cond,1}(:,j,:)] = gaitSPEC_cond.PSC1(:,sub_comp(n),:);

                    [SPEC_subRest{cond,1}(:,j,:)] = gaitSPEC_cond.gaitSPEC_subRest(:,sub_comp(n),:);
                catch 
                    [SPEC_corr{cond,1}(:,j,:)] = nan(1,length(gaitSPEC.freqs));

                    [SPEC_corr_addRest{cond,1}(:,j,:)] = nan(1,length(gaitSPEC.freqs));

                    [SPEC_raw{cond,1}(:,j,:)] = nan(1,length(gaitSPEC.freqs));

                    [PSC1{cond,1}(:,j,:)] = nan(1,length(gaitSPEC.freqs));

                    [SPEC_subRest{cond,1}(:,j,:)] = nan(1,length(gaitSPEC.freqs));

                    fprintf('skipping %s %s \n',sub_ID,all_cond_trials{cond})
                end
            end
            j = j + 1;
        end
        TMP_STUDY.cluster(i).group = group_list(TMP_STUDY.cluster(i).sets);

        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H1000s')) = {1};
        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H2000s')) = {2};
        TMP_STUDY.cluster(i).group(strcmp(TMP_STUDY.cluster(i).group,'H3000s')) = {2};
        TMP_STUDY.cluster(i).group = cellfun(@(x) [[]; x], TMP_STUDY.cluster(i).group);

        group_store = TMP_STUDY.cluster(i).group;

        SPCA_results_spec.SPEC_corr = SPEC_corr;
        SPCA_results_spec.SPEC_corr_addRest = SPEC_corr_addRest;
        SPCA_results_spec.SPEC_raw = SPEC_raw;
        SPCA_results_spec.PSC1 = PSC1;
        SPCA_results_spec.SPEC_subRest = SPEC_subRest;
        SPCA_results_spec.SPEC_rest = SPEC_rest;

        SPCA_results_spec.cluster = i;
        SPCA_results_spec.group_store = group_store;
        SPCA_results_spec.freqs = gaitSPEC.freqs;

        %%%%%%%%%%%%CHANGE CHANGE
        clust_i = CLUSTER_SELECT;
        cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];

        if ~exist(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)]))
            mkdir(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)]))
        end

        save(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)],'SPCA_results_SPEC.mat'),'SPCA_results_spec');


        color.all_trial = [color.terrain;color.speed];
        color.all_trial_shade = [color.terrain_shade;color.speed_shade];
        % - Sanity check figures
        % PSC1
        % figure('color','w');
        % subplot(1,2,1);hold on;
        % for cond = 5:8
        %     x = 0:250;
        %     y = SPCA_results_spec.PSC1{cond,1};
        %     E = std(y,[],2)/sqrt(size(y,2));
        % %     JackKnife_sung(x,mean(y,2),E,color.terrain(cond+1,:),color.terrain_shade(cond+1,:));
        %     plot(x,y,'color',color.all_trial(cond+1,:));
        % end
    end
end
if isunix;exit;end
% keyboard
%% Supplementary figure!! - plot the raw vs. corrected PSDs
% using cluster 7
plotSupplementaryFigure = 1;
if plotSupplementaryFigure 
    
    i = 5;% Posterior Parietal Right
    study_fName = sprintf('temp_study_rejics_spca_%s_%i',study_fName_1, MIN_ICS_SUBJ);
    tmp_dir = [save_dir filesep sprintf('icrej_spca_%s_%i',study_fName_1,MIN_ICS_SUBJ)];
    clust_i = CLUSTER_SWEEP_VALS(3);
    cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)],'SPCA_results_SPEC.mat'));  

    design = 1:4;
    f1 = figure('color','w','position',[100 100 1200 300]);

    subplot(1,5,1)
    for cond = design
        x = 0:250;
        y = SPCA_results_spec.SPEC_raw{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,nanmean(y,2),E,color.all_trial(cond+1,:),color.all_trial(cond+1,:));
    end
    x = 0:250;
    y = SPCA_results_spec.SPEC_rest{1,1};
    E = std(y,[],2)/sqrt(size(y,2));
    JackKnife_sung(x,nanmean(y,2),E,color.all_trial(1,:),color.all_trial_shade(1,:));
    xlim([3 70]);
    xlabel('Frequency (Hz)');
    ylabel('Log Power (dB)');
    title('PSD_{Raw}','color',color.Blue);
    ylim([-35 -10]);
    set(gca,'fontsize',10);
    
    subplot(1,5,2)
    for cond = design
        x = 0:250;
        y = SPCA_results_spec.SPEC_corr_addRest{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,nanmean(y,2),E,color.all_trial(cond+1,:),color.all_trial(cond+1,:));
    end
    x = 0:250;
    y = SPCA_results_spec.SPEC_rest{1,1};
    E = std(y,[],2)/sqrt(size(y,2));
    JackKnife_sung(x,nanmean(y,2),E,color.all_trial(1,:),color.all_trial_shade(1,:));
    xlim([3 70]);ylim([-35 -10]);
    xlabel('Frequency (Hz)');
    ylabel('');
    title('PSD_{sPCA corrected} + PSD_{Rest}','color',color.Red);
    set(gca,'fontsize',10);
    
    subplot(1,5,3)
    for cond = design
        x = 0:250;
        y = SPCA_results_spec.SPEC_subRest{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,nanmean(y,2),E,color.all_trial(cond+1,:),color.all_trial(cond+1,:));
    end
    xlim([3 70]);tmp_gca = gca; ylim_2 = tmp_gca.YLim;
    xlabel('Frequency (Hz)');
    ylabel('');
    title('PSD_{Raw} - PSD_{Rest}','color',color.Green);
    set(gca,'fontsize',10);
    
    subplot(1,5,4)
    for cond = design
        x = 0:250;
        y = SPCA_results_spec.SPEC_corr{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,nanmean(y,2),E,color.all_trial(cond+1,:),color.all_trial(cond+1,:));
    end
    xlim([3 70]);ylim(ylim_2);
    xlabel('Frequency (Hz)');
    ylabel('');
    title('PSD_{sPCA corrected}','color',color.Yellow);
    set(gca,'fontsize',10);
    
    subplot(1,5,5)
    for cond = design
        x = 0:250;
        y = SPCA_results_spec.PSC1{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,nanmean(y,2),E,color.all_trial(cond+1,:),color.all_trial(cond+1,:));
    end
    xlim([3 70]);
    xlabel('Frequency (Hz)');
    ylabel('');
    title('PSC removed');
    set(gca,'fontsize',10);
%     saveas(f1,fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(i)],'SPCA_results_SPEC.jpg'));  
    %}

    %% Sanity check - compare corrected vs. original PSD
    figure('color','w');
    subplot(2,2,1)
    for cond = 1:4
        x = 0:250;
        y = SPCA_results_spec.SPEC_corr_addRest{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,mean(y,2),E,color.all_trial(cond+1,:),color.all_trial_shade(cond+1,:));
    end
    xlim([3 70]);ylim([-36 -10])
    xlabel('Frequency (Hz)');
    ylabel('Log Power (dB)');
    title('SPEC Corrected');

    subplot(2,2,2)
    for cond = 1:4
        x = 0:250;
        y = SPCA_results_spec.SPEC_raw{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,mean(y,2),E,color.all_trial(cond+1,:),color.all_trial_shade(cond+1,:));
    end
    xlim([3 70]);ylim([-36 -10])
    xlabel('Frequency (Hz)');
    ylabel('Log Power (dB)');
    title('SPEC Original');

    subplot(2,2,3)
    for cond = 5:8
        x = 0:250;
        y = SPCA_results_spec.SPEC_corr_addRest{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,mean(y,2),E,color.all_trial(cond+1,:),color.all_trial_shade(cond+1,:));
    end
    xlim([3 70]);ylim([-36 -10])
    xlabel('Frequency (Hz)');
    ylabel('Log Power (dB)');
    title('SPEC Corrected');

    subplot(2,2,4)
    for cond = 5:8
        x = 0:250;
        y = SPCA_results_spec.SPEC_raw{cond,1};
        E = std(y,[],2)/sqrt(size(y,2));
        JackKnife_sung(x,mean(y,2),E,color.all_trial(cond+1,:),color.all_trial_shade(cond+1,:));
    end
    xlim([3 70]);ylim([-36 -10])
    xlabel('Frequency (Hz)');
    ylabel('Log Power (dB)');
    title('SPEC Original');
end
%% Gather Raw SPEC PSD (not SPCA corrected) for each cluster
% Cautious: take a long time to run, it is not used anymore
% if ~ispc
%     cluster_dir = convertPath2UNIX(CLUSTER_DIRS{k_i});
% else
%     cluster_dir = convertPath2Drive(CLUSTER_DIRS{k_i});
% end
generateOriginalPSD = false
if generateOriginalPSD
    plot_store_dir = [cluster_dir filesep 'plots_out'];
    if ~exist(plot_store_dir,'dir')
        mkdir(plot_store_dir);
    end
    spec_data_dir = [cluster_dir filesep 'in_brain' filesep 'spec_data'];
    if ~exist(spec_data_dir,'dir')
        mkdir(spec_data_dir)
    end

    mim_gen_cluster_spec_CL(TMP_STUDY,TMP_ALLEEG,cluster_dir,spec_data_dir,CLUSTER_PICKS);
    mim_gen_cluster_spec_CL(TMP_STUDY_gait,TMP_ALLEEG_gait,cluster_dir,spec_data_dir,CLUSTER_PICKS);
end
%% =================== 
%% FOOOF Analysis
% j = 4;
% reduce_method = 'max_iclabel';
% CLUSTER_SELECT = CLUSTER_SWEEP_VALS(j);
% clust_i = CLUSTER_SELECT;
% cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];
% load([cluster_dir filesep 'ERSP_Plots' filesep 'cluster_info.mat'])
%% Generate spectral 
SUB_GROUP_FNAME_REGEX = [];
% Setup Fooof
SAVE_DATA = true;
settings = struct();  % Use defaults
settings.peak_width_limits = [1 8];%default [1 6] % Amanda used [1 8] in her paper
settings.min_peak_height = 0.05;
% settings.peak_threshold = 2;
settings.max_n_peaks = 3; % originally set to be 2 - 2023-06-07
% the settings are consitent with fooof on github
f_range = [3, 40];
theta_band = [4, 8];
alpha_band = [8 12];
beta_band  = [12 30];
cond_terrains = {'flat','low','med','high','rest'};
cond_speeds = {'0p25','0p5','0p75','1p0','rest'};
COND_CHARS = {cond_terrains,cond_speeds};
% spec_data_cl2_0p250p50p751p0_subbaselined_commonbase
% spec_data_cl2_flatlowmedhigh_subtractmean
keywords = {'Terrain'};

spec_data_dir = [cluster_dir filesep 'in_brain' filesep 'spec_data'];
save_spec_dir = [spec_data_dir filesep 'psd_calcs_spca'];
if ~exist(save_spec_dir,'dir')
    mkdir(save_spec_dir);
end

%% READ IN SUBJECT SPECIFIC SPEEDS FOR TERRAIN
SPEED_CUTOFF = 0.1;
MasterTable = mim_read_master_sheet();
speed_table = table(categorical(MasterTable.subject_code),MasterTable.terrain_trials_speed_ms);
speed_alleeg = cell(length(ALLEEG),2);
for i = 1:length(ALLEEG)
    ss = ALLEEG(i).subject;
    ind = speed_table.Var1==ss;
    chk1 = strcmp(ALLEEG(i).group,SUB_GROUP_FNAME_REGEX) || isempty(SUB_GROUP_FNAME_REGEX);
    if any(ind) && chk1
        speed_alleeg{i,1} = speed_table.Var1(ind);
        speed_alleeg{i,2} = double(speed_table.Var2(ind));
    end
end
speed_alleeg = speed_alleeg(~cellfun(@isempty,speed_alleeg(:,1)),:);

%% (TABLE) GENERATE FOOOF VALUES ======================================= %%
cluster_dir = [tmp_dir filesep num2str(CLUSTER_SELECT) filesep reduce_method];
DESIGN_I = 1;
fooof_group_results_org = cell(1,length(DESIGN_I));
fooof_results = cell(length(DESIGN_I),1);
fooof_diff_store = cell(length(DESIGN_I),1);
fooof_apfit_store = cell(length(DESIGN_I),1);
spec_data_original = cell(length(DESIGN_I),1);

for g = DESIGN_I
    for k = valid_clusters
        fprintf('Processing %s \n',num2str(k));
        cluster_update = load([cluster_dir filesep 'in_brain' filesep sprintf('cl_inf_spca_inbrain_%s_%i_%s.mat',study_fName_1,CLUSTER_SELECT,reduce_method)]); %par_load(cluster_dir,sprintf('cluster_inf_%i.mat',clust_i));
        cluster_update = cluster_update.SAVEVAR;
        TMP_STUDY.cluster = cluster_update;
        
        file_mat = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results_SPEC.mat');        
        tmp_spec = load(file_mat);

        % Note: spec_subj_mean_stats separate young and older adults
        specdata_rest = tmp_spec.SPCA_results_spec.SPEC_rest';
        specfreqs_rest = tmp_spec.SPCA_results_spec.freqs';

        % spca corrected
        specdata_gait_spca = tmp_spec.SPCA_results_spec.SPEC_corr_addRest';
        specfreqs_gait_spca = tmp_spec.SPCA_results_spec.freqs';

        % original psd
        specdata_gait_raw = tmp_spec.SPCA_results_spec.SPEC_raw';
        specfreqs_gait_raw = tmp_spec.SPCA_results_spec.freqs';

        if g == 1 %design = 1 terrain
            [specdata_all_2group_spca,specdata_all_3group_spca] = organize_specdata_2(specdata_rest,specdata_gait_spca,TMP_STUDY,k,1:4);
            [specdata_all_2group_raw,specdata_all_3group_raw] = organize_specdata_2(specdata_rest,specdata_gait_raw,TMP_STUDY,k,1:4);
        elseif g == 2
            [specdata_all_2group_spca,specdata_all_3group_spca] = organize_specdata_2(specdata_rest,specdata_gait_spca,TMP_STUDY,k,5:8);
            [specdata_all_2group_raw,specdata_all_3group_raw] = organize_specdata_2(specdata_rest,specdata_gait_raw,TMP_STUDY,k,5:8);
        end

        clear cfg;
        cfg.f_range = f_range;
        cfg.settings = settings;
        cfg.specfreqs_rest = specfreqs_rest;
        fooof_outcome_spca = fooof_process(specdata_all_2group_spca, TMP_STUDY,cluster_update,k,g,cfg);
        fooof_group_results_org_spca{g}{k} = fooof_outcome_spca.fooof_group_results_org{g}{k};
        fooof_diff_store_spca{g}{k} = fooof_outcome_spca.fooof_diff_store{g}{k};
        fooof_apfit_store_spca{g}{k} = fooof_outcome_spca.fooof_apfit_store{g}{k};
        spec_data_original_spca{g}{k} = fooof_outcome_spca.spec_data_original{g}{k};

        % run fooof on non-spca corrected
        fooof_outcome_raw = fooof_process(specdata_all_2group_raw, TMP_STUDY,cluster_update,k,g,cfg);
        fooof_group_results_org_raw{g}{k} = fooof_outcome_raw.fooof_group_results_org{g}{k};
        fooof_diff_store_raw{g}{k} = fooof_outcome_raw.fooof_diff_store{g}{k};
        fooof_apfit_store_raw{g}{k} = fooof_outcome_raw.fooof_apfit_store{g}{k};
        spec_data_original_raw{g}{k} = fooof_outcome_raw.spec_data_original{g}{k};
    
        % fooof ID
        num_young = size(specdata_all_2group_spca{1},2); 
        num_old = size(specdata_all_2group_spca{2},2);
        fooof_subID{g}{k}{1} = cluster_update(k).sets(1:num_young);
        fooof_subID{g}{k}{2} = cluster_update(k).sets(num_young+1:end);
        
        cl_chars = {TMP_STUDY.datasetinfo(cluster_update(k).sets).subject};
        fooof_subj_char{g}{k}{1} = {cl_chars{1:num_young}};
        fooof_subj_char{g}{k}{2} = {cl_chars{num_young+1:end}};
    end
end

% add subjectID 


save([save_spec_dir filesep 'fooof_results_20250511.mat'],'fooof_outcome_spca','fooof_group_results_org_spca',...
    'fooof_diff_store_spca','fooof_apfit_store_spca','spec_data_original_spca',...
    'fooof_group_results_org_raw','fooof_diff_store_raw','fooof_apfit_store_raw',...
    'spec_data_original_raw','fooof_subj_char','fooof_subID');

% organize
save_filename.xlsx_name = 'fooof_group_results_org_spca_20250511';
save_filename.psd_table_name = 'psd_feature_table_spca_20250511';
save_filename.psd_table_unstack_name = 'psd_feature_table_unstack_spca_20250511';
save_filename.save_dir = save_spec_dir;
[psd_feature_table_spca,psd_feature_table_unstack_spca,psd_feature_table_simplify_spca] = fooof_result_organize(fooof_group_results_org_spca,DESIGN_I,save_filename);


