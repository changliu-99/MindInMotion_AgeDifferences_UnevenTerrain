% Calculate peak-to-peak ERSP, and perform statistical analsysis of ERSP results
% Chang Liu - 20240614


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
cd(run_dir)
%- addpath for local folder
addpath(source_dir)
addpath(run_dir)
addpath(PATH_ROOT)
addpath([filesep 'blue' filesep 'dferris' filesep USER_NAME filesep 'scripts' filesep 'MiM_CRUNCH' filesep '_submodules' filesep 'hline_vline']);
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
CLUSTER_SWEEP_VALS = 11;%[12,14,19]; %[10,13,14,19,20]; %K_RANGE(1):K_RANGE(2);
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

%% Setup STUDY stats
%* ERSP PARAMS
STUDY = [];
STUDY.etc = [];

ERSP_STAT_PARAMS = struct('condstats','on',... % ['on'|'off]
    'groupstats','off',... %['on'|'off']
    'method','bootstrap',... % ['param'|'perm'|'bootstrap']
    'singletrials','off',... % ['on'|'off'] load single trials spectral data (if available). Default is 'off'.
    'mode','fieldtrip',... % ['eeglab'|'fieldtrip']
    'fieldtripalpha',0.05,... % [NaN|alpha], Significance threshold (0<alpha<<1)
    'fieldtripmethod','montecarlo',... %[('montecarlo'/'permutation')|'parametric']
    'fieldtripmcorrect','cluster',...  % ['cluster'|'fdr']
    'fieldtripnaccu',4000);

STUDY = pop_statparams(STUDY,'condstats',ERSP_STAT_PARAMS.condstats,...
        'groupstats',ERSP_STAT_PARAMS.groupstats,...
        'method',ERSP_STAT_PARAMS.method,...
        'singletrials',ERSP_STAT_PARAMS.singletrials,'mode',ERSP_STAT_PARAMS.mode,...
        'fieldtripalpha',ERSP_STAT_PARAMS.fieldtripalpha,...
        'fieldtripmethod',ERSP_STAT_PARAMS.fieldtripmethod,...
        'fieldtripmcorrect',ERSP_STAT_PARAMS.fieldtripmcorrect,'fieldtripnaccu',ERSP_STAT_PARAMS.fieldtripnaccu);

ERSP_PARAMS = struct('subbaseline','off',...
    'timerange',[],...
    'ersplim',[-2,2],...
    'freqfac',4,...
    'cycles',[3,0.8],...
    'freqrange',[1,200],...
    'plot_freqrange',[4,60],...
    'plot_clim',[-2,2]);
STUDY.etc.erspparams = ERSP_PARAMS;

study_fName = sprintf('temp_study_rejics_spca_%s_%i',study_fName_1, MIN_ICS_SUBJ);
tmp_dir = [save_dir filesep sprintf('icrej_spca_%s_%i',study_fName_1,MIN_ICS_SUBJ)];

TMP_STUDY = importdata([tmp_dir filesep sprintf('temp_study_rejics_spca_%s_%i.study',study_fName_1_all, MIN_ICS_SUBJ)]);
TMP_STUDY.design = [];
TMP_STUDY = std_makedesign(TMP_STUDY, [], 1, 'subjselect',{},'variable1','cond','values1', {'flat','low','med','high'});
TMP_STUDY.design.variable(2) = [];
STUDY.design = TMP_STUDY.design;
STUDY.currentdesign = 1;

TMP_STUDY_group = std_makedesign(TMP_STUDY, [], 1, 'subjselect',{},'variable1','cond','values1', {'flat','low','med','high'},...
    'variable2','group');
TMP_STUDY_group.design.variable(2).value = {'1','2'};
%% Setup ersp comparison
fcn = @erspStats;
fcn2 = @erspStats_V2;

%% ===================================================================== %%
%% Load SPCA results
CLUSTER_SELECT = 11;
reduce_method = 'max_iclabel';
clust_i = CLUSTER_SELECT;
cluster_dir = [tmp_dir filesep num2str(clust_i) filesep reduce_method];
cluster_update = load([cluster_dir filesep 'in_brain' filesep sprintf('cl_inf_spca_inbrain_%s_%i_%s.mat',study_fName_1,clust_i,reduce_method)]); %par_load(cluster_dir,sprintf('cluster_inf_%i.mat',clust_i));
TMP_STUDY.cluster = cluster_update.SAVEVAR;

[~,main_cl_inds,~,valid_clusters,~,nonzero_clusters] = eeglab_get_cluster_comps_CL(TMP_STUDY);

load(fullfile(load_dir,'H1004','GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh','H1004.icatimef'),'-mat','parameters');
warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
group = 2;
%%
% old reduction method 
% if CLUSTER_SELECT == 13
%     cluster_label = {'SuppMotor','Cingulate','SMA','Occipital','Caudate','SMA','Cingulate','SuppMotor','PP','PP'};
%     climMax = [0.5 0.5 1.2 0.5 0.5 1.2 0.5 0.5 1 1];
% elseif CLUSTER_SELECT == 12
%     cluster_label = {'SMA','Cingulate','PP','Occipital','SuppMotor','SMA','Cingulate','Cingulate','PP','Occipital'};
%     climMax = [1 0.5 0.5 0.5 0.5 1 0.5 0.5 0.5 1];
% elseif CLUSTER_SELECT == 11
%     cluster_label = {'Occipital','','SuppMotor','Cingulate','PP','SuppMotor','PP','SMA','SMA'};
%     climMax = [0.5 0.5 0.5 0.5 1 0.5 1 1 1];
% end

switch reduce_method
    case 'max_iclabel'
        if CLUSTER_SELECT == 11
            cluster_name = {'SMA','SuppMotor','PP','SuppMotor','PP','SMA','PP','SuppMotor',''};
            climMax = [1 0.5 0.5 0.5 0.5 1 0.5 0.5 0.5];
        end
end

config.warpingvalues = warpingvalues;
config.colormap_ersp = colormap_ersp;
config.title_keyword = {'flat','low','med','high'};
config.STUDY = STUDY;
config.runPairwise = 0;
config.YlimMax = 50;
config.CLUSTER_SWEEP_VALS = CLUSTER_SELECT;
config.color = color;

p2p_ERSP_table_sum = table;

for k = valid_clusters
   
    config.k = k;
    %% LOAD SPCA result
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.mat'));
    sub_ID = TMP_STUDY.subject{1};
    load_ersp_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitERSP_',sub_ID,'_','flat']);
    load(load_ersp_spca);
    SPCA_results.times = gaitERSP_cond.times;
    SPCA_results.logfreqs = gaitERSP_cond.logfreqs;
    
    % First step, run stats for the GPM style data
    allersp = SPCA_results.ERSP_GPM_corr;
    allersp_ERSP = SPCA_results.ERSP_corr;
    for j = 1:4
        allersp_reorder{j,1} = permute(allersp{j},[3 1 2]);%transform to freq x time x subject
    end
    % organize the allersp to young and older adults
    for j = 1:4
        allersp_reorder_group{j,1} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 1);
        allersp_reorder_group{j,2} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 2);
    end
    % - Compute the peak to peak ERSP
    config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'p2p-ERSP_group.jpg');
    config.save_file = 0;
    config.climMat_max = climMax(find(valid_clusters == k));
    plot_compute_PtP_ERSP(allersp_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
    
    config.save_file = 0;
    p2p_ERSP_table = compute_PtP_ERSP(allersp_reorder,SPCA_results.times,SPCA_results.logfreqs,config);
    p2p_ERSP_table = compute_PtP_ERSP(allersp_reorder,SPCA_results.times,SPCA_results.logfreqs,config);
    
    subj = TMP_STUDY.subject;
    p2p_ERSP_table_all = [table(repmat(TMP_STUDY.cluster(k).sets',16,1),'VariableNames',{'subID'})... 
            table(repmat(subj(TMP_STUDY.cluster(k).sets)',16,1),'VariableNames',{'subjectName'}) ...
            table(repmat(k,height(p2p_ERSP_table),1),'VariableNames',{'cluster'}) ...
            table(repmat(SPCA_results.group_store',16,1),'VariableNames',{'group'}) p2p_ERSP_table];

    
    p2p_ERSP_table_sum = [p2p_ERSP_table_sum;p2p_ERSP_table_all];
            
end
saveP2Ptable = 0;
if saveP2Ptable
    p2p_ERSP_table_sum_unstack = unstack(p2p_ERSP_table_sum,{'p2p_ersp'},'band','VariableNamingRule','preserve');

    % -- SAVE simplified table that only contained relevant variables - update 20240425
    spec_data_dir = [cluster_dir filesep 'in_brain' filesep 'spec_data'];
    save_dir = [spec_data_dir filesep 'psd_calcs_spca'];
    mkdir(save_dir)

    save_filename.xlsx_name = 'p2p_ERSP_table_sum_202505';
    save_filename.psd_table_name = 'p2p_ERSP_table_sum_202505';
    save_filename.save_dir = save_dir;

    writetable(p2p_ERSP_table_sum,[save_filename.save_dir filesep [save_filename.xlsx_name,'.xlsx']]);
    save([save_filename.save_dir filesep [save_filename.psd_table_name,'.mat']],'p2p_ERSP_table_sum');

    writetable(p2p_ERSP_table_sum_unstack,[save_filename.save_dir filesep [save_filename.xlsx_name,'_unstack.xlsx']]);
    save([save_filename.save_dir filesep [save_filename.psd_table_name,'_unstack.mat']],'p2p_ERSP_table_sum');
end

%% Paper Figure 7-9: Peak-to-peak update: corrected mean for each participant and then run peak-to-peak again
% did both median and mean across participants. did not see much
% differences in the plot
config.warpingvalues = warpingvalues;
config.colormap_ersp = colormap_ersp;
config.title_keyword = {'flat','low','med','high'};
config.STUDY = STUDY;
config.runPairwise = 0;
config.YlimMax = 50;
config.CLUSTER_SWEEP_VALS = CLUSTER_SELECT;
config.color = color;
p2p_ERSP_table_correctMean_sum = table;
swing_ERSP_table_all_sum = table;
for k = valid_clusters
    config.k = k;
    %% LOAD SPCA result
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.mat'));
    sub_ID = TMP_STUDY.subject{1};
    load_ersp_spca = fullfile(load_dir, sub_ID,'GAIT_EPOCHED_ALL','0p250p50p751p0flatlowmedhigh',['gaitERSP_',sub_ID,'_','flat']);
    load(load_ersp_spca);
    SPCA_results.times = gaitERSP_cond.times;
    SPCA_results.logfreqs = gaitERSP_cond.logfreqs;
    
    % First step, run stats for the GPM style data
    allersp = SPCA_results.ERSP_GPM_corr;
    allersp_ERSP = SPCA_results.ERSP_corr;
    for j = 1:4
        allersp_reorder{j,1} = permute(allersp{j},[3 1 2]);%transform to freq x time x subject
    end
    % organize the allersp to young and older adults
    for j = 1:4
        allersp_reorder_group{j,1} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 1);
        allersp_reorder_group{j,2} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 2);
    end
    config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'p2p-ERSP_group_correctMean.tif');
    config.save_pdf_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'p2p-ERSP_group_correctMean.pdf');
    config.save_file = 1;
    config.climMat_max = climMax(find(valid_clusters == k));
    config.cond = 1;
    p2p_ERSP_table_correctMean = plot_compute_PtP_ERSP_correctMean(allersp_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
    
    subj = TMP_STUDY.subject;
    p2p_ERSP_table_correctMean_all = [table(repmat(TMP_STUDY.cluster(k).sets',16,1),'VariableNames',{'subID'})... %3 band 4 conditions
            table(repmat(subj(TMP_STUDY.cluster(k).sets)',16,1),'VariableNames',{'subjectName'}) ...
            table(repmat(k,height(p2p_ERSP_table_correctMean),1),'VariableNames',{'cluster'}) ...
            table(repmat(SPCA_results.group_store',16,1),'VariableNames',{'group'}) p2p_ERSP_table_correctMean];

    p2p_ERSP_table_correctMean_sum = [p2p_ERSP_table_correctMean_sum;p2p_ERSP_table_correctMean_all];
   
   %% Evaluate alpha and beta desynchronization during contralateral swing 
    swing_ERSP_table = plot_compute_contralateral_swing_ersp(allersp_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
    subj = TMP_STUDY.subject;
    swing_ERSP_table_all = [table(repmat(TMP_STUDY.cluster(k).sets',16,1),'VariableNames',{'subID'})... %3 band 4 conditions
            table(repmat(subj(TMP_STUDY.cluster(k).sets)',16,1),'VariableNames',{'subjectName'}) ...
            table(repmat(k,height(swing_ERSP_table),1),'VariableNames',{'cluster'}) ...
            table(repmat(SPCA_results.group_store',16,1),'VariableNames',{'group'}) swing_ERSP_table];
    swing_ERSP_table_all_sum = [swing_ERSP_table_all_sum;swing_ERSP_table_all];
end
save_p2p_correctMean_table = 0; %SAVE the p2p ersp data
if save_p2p_correctMean_table
    % -- SAVE the p2p table
    save_data_dir = ['M:\liu.chang1\scripts\MiM_CRUNCH\5_stats'];

    save_filename.xlsx_name = 'p2p_correctMean_table_2025_05';
    save_filename.psd_table_name = 'p2p_correctMean_table_2025_05';
    save_filename.save_dir = save_data_dir;

    writetable( p2p_ERSP_table_correctMean_sum,[save_filename.save_dir filesep [save_filename.xlsx_name,'.xlsx']]);
    save([save_filename.save_dir filesep [save_filename.psd_table_name,'.mat']],'p2p_ERSP_table_correctMean_sum');
end
keyboard
%% Paper figure: Violin plot of peak to peak ERSP
data_table = p2p_ERSP_table_correctMean_sum;
outcome_measure = {'p2p_ersp'};
% data_table = swing_ERSP_table_all_sum;
% outcome_measure = 'LTO_LHS_ERSP';% - right sensorimotor 13
% outcome_measure = 'LTO_RTO_ERSP';% - right sensorimotor 13
% outcome_measure = 'RTO_RHS_ERSP';% - left sensorimotor 11
% outcome_measure = 'RTO_LTO_ERSP';% - left sensorimotor 11

band = {'\theta','\alpha','\beta'};
title_plot = {'\theta','\alpha','\beta'}
color_dark = color.terrain(2:end,:);
color_light = color.terrain_old(2:end,:);
xtick_label_g = {'rest','flat','low','med','high'};
ylabel_name = {''};


for k = valid_clusters
    p = 1;
    figure('color','white','unit','centimeters','position',[10 10 6.5 20]);
    for j = 1:length(band)
        for i = 1
            fig = subplot(3,1,p);
            hold on;
            data_plot = data_table(data_table.cluster == k & strcmp(data_table.band,band{j}),:);
            % there is a huge outlier for k = 8 at theta band,remove that
            % person
%             if k == 8
%                 idx1 = find(data_plot.p2p_ersp > 5.52); %mean+-5std;
%                 data_plot.p2p_ersp(idx1) = NaN;
%                 data_plot.min_ersp(idx1) = NaN;
%                 data_plot.max_ersp(idx1) = NaN;
%             end
            
            hold on;
            f1 = violinplot(data_plot.(outcome_measure{i})(data_plot.group == 1),double(data_plot.cond(data_plot.group == 1)),...
                'ViolinColor',{color_dark(:,:)},'ViolinAlpha',{0.2 0.3},...
                'MarkerSize',10,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter',...
                'width',0.3,'HalfViolin','left',...
                'ShowMean', false,'BoxColor',[0 0 0]);
    %         for cond = 1:4
    %             mean_f1(cond) = median(data_plot.p2p_ersp(data_plot.group == 1 & data_plot.cond == cond));
    %             mean_f2(cond) = median(data_plot.p2p_ersp(data_plot.group == 2 & data_plot.cond == cond));
    %         end
            dataObjs = findobj( fig,'-property','YData');
            for ii = 1:length(dataObjs)
                dataObjs(ii).XData = dataObjs(ii).XData -0.2;
            end
            f2 = violinplot(data_plot.(outcome_measure{i})(data_plot.group == 2),double(double(data_plot.cond(data_plot.group == 2)))+1.5,...
                'ViolinColor',{color_light(:,:)},'ViolinAlpha',{0 0},...
                'MarkerSize',10,'EdgeColor',[0.5 0.5 0.5],'DataStyle', 'scatter',...
                'width',0.3,'HalfViolin','right',...
                'ShowMean', false,'BoxColor',[0 0 0]);
    %             ylim([0 3]);
            if i == 1;ylabel(sprintf('%s Power(dB)',title_plot{j}),'fontweight','bold');end

    %         plot(unique(double(data_plot.cond(data_plot.group == 1)))-0.2,mean_f1,'-','color',[174,1,126]/255,'linewidth',1.5);
    %         plot(unique(double(data_plot.cond(data_plot.group == 2))),mean_f2,'-','color',[8,81,156]/255,'linewidth',1.5);
            % axis padded
            fig_i = get(groot,'CurrentFigure');
            % fig_i.Position = [200,200,1820,920];
            box off
            if j == 1;title(ylabel_name{i});end
            if j == 3
                set(gca,'xticklabel', xtick_label_g(2:end),'fontsize',12);
            else
                set(gca,'xticklabel',{'','',''},'fontsize',12);
            end
            xtickangle(45)            
            
            p = p + 1;
        end
                % perform stats
        data_plot.cond = data_plot.cond - 2;
        data_plot.group = data_plot.group - 1;
        data_plot.group = categorical(data_plot.group);
        data_plot.cond = categorical(data_plot.cond);
        
%         data_plot.group = double(data_plot.group);
%         data_plot.cond = double(data_plot.cond);

%!!!!!%%%%%%%%%%%%%%%%%%%
% ANOVA result with interaction term does not match with R !! USE R
% instead!
%%%%%%%%%%%%%%%%%%%%%%%%%
%         data_plot.subID = categorical(data_plot.subID);
%         lm_terrain_p2p = fitlme(data_plot,'p2p_ersp ~ 1+group+cond+group:cond+(1|subjectName)', 'FitMethod','REML');
%         anova(lm_terrain_p2p,'dfmethod','satterthwaite')
%         keyboard
    end
    savedir = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)]);
    exportgraphics(gcf,fullfile(savedir,['cluster_',num2str(k),'_','p2p_ERSP_violin','.tif']))
%     exportgraphics(gcf,fullfile(savedir,['cluster_',num2str(k),'_','p2p_ERSP_violin','.pdf']))
end

%}
%% compute all stats related with ERSPs
sub_ID = TMP_STUDY.subject{1};
for k = valid_clusters
    fprintf('Processing cluster = %s \n',num2str(k));
    performCorrect = true;
    config.k = k;
    %% LOAD SPCA result
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'SPCA_results.mat'));
    sub_ID = TMP_STUDY.subject{1};
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
    
    % organize the allersp to young and older adults
    for j = 1:4
        allersp_reorder_group{j,1} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 1);
        allersp_reorder_group{j,2} = allersp_reorder{j,1}(:,:,SPCA_results.group_store == 2);
        
        allersp_ERSP_reorder_group{j,1} = allersp_ERSP_reorder{j,1}(:,:,SPCA_results.group_store == 1);
        allersp_ERSP_reorder_group{j,2} = allersp_ERSP_reorder{j,1}(:,:,SPCA_results.group_store == 2);
    end
    %% -
%     keyboard
    performCorrect_within = 1;
    if performCorrect_within                
        % Baseline correction within the gait cycle
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_within_young_202505.tif');
        config.save_pdf_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_within_young_202505.pdf');
        config.YlimMax = 50;
        config.climMat_max = climMax(find(valid_clusters == k));
        PerformBaselineCorrect_within(allersp_ERSP_reorder_group(1:4,1),SPCA_results.times,SPCA_results.logfreqs,config)
        exportgraphics(gcf,config.save_filepath);
        exportgraphics(gcf,config.save_pdf_filepath);
        
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_within_old_202505.tif');
        config.save_pdf_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_within_old_202505.pdf');
        config.YlimMax = 50;
        config.climMat_max = climMax(find(valid_clusters == k));
        PerformBaselineCorrect_within(allersp_reorder_group(1:4,2),SPCA_results.times,SPCA_results.logfreqs,config)
        exportgraphics(gcf,config.save_filepath);
        exportgraphics(gcf,config.save_pdf_filepath);
    end         
%         % Common baseline subtraction
%         % How's common baseline subtraction work? Does it subtract the average
%         % for each person or just all
%         % https://sccn.ucsd.edu/pipermail/eeglablist/2012/004628.html
%         common baseline issue
%         config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_common_young.jpg');
%         PerformBaselineCorrect_common(allersp_ERSP_reorder_group(1:4,1),SPCA_results.times,SPCA_results.logfreqs,config)
% 
%         config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],'allersp_spca_common_old.jpg');
%         PerformBaselineCorrect_common(allersp_ERSP_reorder_group(1:4,2),SPCA_results.times,SPCA_results.logfreqs,config)

        
        %% -- Common baseline subtract - Cond only
    performCorrect_common = 1;
    if performCorrect_common        
        % old reduction method 
%         
        % NEW reduction method - max_iclabel
        
        figure_keyword = {'YA','OA'};
        STUDY_V2 = STUDY;
        STUDY_V2.design = TMP_STUDY_group.design;
        STUDY_V2.etc.statistics.groupstats = 'off';
        STUDY_V2.design.variable(2) = [];
        config.STUDY = STUDY_V2;
        config.computeGroupStats = false;
        config.runPairwiseGroup = false;
        config.runPairwiseCond = false;
        config.plotGroup = 1;
        config.climMat_max = climMax(find(valid_clusters == k));
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_',figure_keyword{config.plotGroup},'.jpg']);
        [allersp_out_YA,pcond_ersp_crop_YA, ~, ~] = PerformBaselineCorrect_common_CondOnly(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);

        
        figure_keyword = {'YA','OA'};
        STUDY_V2 = STUDY;
        STUDY_V2.design = TMP_STUDY_group.design;
        STUDY_V2.etc.statistics.groupstats = 'off';
        STUDY_V2.design.variable(2) = [];
        config.STUDY = STUDY_V2;
        config.computeGroupStats = false;
        config.runPairwiseGroup = false;
        config.runPairwiseCond = false;
        config.plotGroup = 2;
        config.climMat_max = climMax(find(valid_clusters == k));
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_',figure_keyword{config.plotGroup},'.jpg']);
        [allersp_out_OA,pcond_ersp_crop_OA, ~, ~] = PerformBaselineCorrect_common_CondOnly(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);

        %% --- Common baseline subtract - Group and Cond together
        % fix of 2way ANOVA: https://sccn.ucsd.edu/pipermail/eeglablist/2022/016748.html
        % it is not recommended to get interaction effect
        figure_keyword = {'YA','OA'};
        STUDY_V2 = STUDY;
        STUDY_V2.design = TMP_STUDY_group.design;
        STUDY_V2.etc.statistics.groupstats = 'on';
        
        config.STUDY = STUDY_V2;
        config.computeGroupStats = true;
        config.runPairwiseGroup = true;
        config.runPairwiseCond = false;
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.jpg']);
        [allersp_out_YAOA,pcond_ersp_crop_YAOA, pgroup_ersp_crop_YAOA, pinter_ersp_crop_YAOA] = PerformBaselineCorrect_common(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);

        config.STUDY = STUDY_V2;
        config.computeGroupStats = true;
        config.runPairwiseGroup = true;
        config.runPairwiseCond = false;
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_pairwise','.jpg']);
        runPairwiseGroup(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
    end
        %% ERSP relative to flat
    performCorrect_rel_flat = 1;
    if performCorrect_rel_flat   
        % first make the data relative to flat
        STUDY_V2 = STUDY;
        STUDY_V2.design = TMP_STUDY_group.design;
        STUDY_V2.etc.statistics.groupstats = 'on';

        config.STUDY = STUDY_V2;
        config.computeGroupStats = true;
        config.runPairwiseGroup = true;
        config.runPairwiseCond = true;
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_pairwise_rel_flat','.pdf']);
        [erspDiff_wind] = runPairwiseCond_rel_flat(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
        save(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_pairwise_rel_flat.mat']),'erspDiff_wind');
        
        config.runPairwiseGroup = true;
        config.runPairwiseCond = false;
        config.save_filepath = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_pairwise_group_rel_flat','.pdf']);
        [erspDiff_wind_group] = runPairwiseGroup_rel_flat(allersp_ERSP_reorder_group,SPCA_results.times,SPCA_results.logfreqs,config);
        save(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_pairwise_group_rel_flat.mat']),'erspDiff_wind');

    end

end
% what about common baseline and CRUNCH
if isunix; exit;end;
%%  Plot only 
close all
% compute all stats
switch reduce_method
    case 'max_iclabel'
        if CLUSTER_SELECT == 11
            cluster_name = {'SMA','SuppMotor','PP','SuppMotor','PP','SMA','PP','SuppMotor',''};
            climMax = [1 0.5 1 0.5 1 1 1 0.5 0.5];
        end
end

for k = valid_clusters(1)
    config.k = k;
    config.climMat_max = climMax(find(valid_clusters == k));
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.mat']));  
    config.save_filepath =  fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.jpg']);
    FIGURE_performBaselineCorrect_common(allerspdata_crop,SPCA_results.times,SPCA_results.logfreqs,pcond_ersp_crop,pgroup_ersp_crop,pinter_ersp_crop,config)
    
%     mkdir(fullfile('J:\ChangLiu\UFL Dropbox\Chang Liu\0.0 Writing_postdoc\MiM_CRUNCH\Figures',['Cluster_',num2str(k)]));
%     config.save_filepath_dropbox = fullfile('J:\ChangLiu\UFL Dropbox\Chang Liu\0.0 Writing_postdoc\MiM_CRUNCH\Figures',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.pdf']);
%     exportgraphics(gcf,config.save_filepath_dropbox);
    config.save_filepath_M = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.pdf']);
    exportgraphics(gcf,config.save_filepath_M);
end

for k = valid_clusters(1)
    config.k = k;
    config.climMat_max = climMax(find(valid_clusters == k));
    load(fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond','.mat']));  
    config.save_filepath =  fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond_tworow','.jpg']);
    FIGURE_performBaselineCorrect_common_tworow(allerspdata_crop,SPCA_results.times,SPCA_results.logfreqs,pcond_ersp_crop,pgroup_ersp_crop,pinter_ersp_crop,config)
    
%     mkdir(fullfile('J:\ChangLiu\UFL Dropbox\Chang Liu\0.0 Writing_postdoc\MiM_CRUNCH\Figures',['Cluster_',num2str(k)]));
%     config.save_filepath_dropbox = fullfile('J:\ChangLiu\UFL Dropbox\Chang Liu\0.0 Writing_postdoc\MiM_CRUNCH\Figures',['Cluster_',num2str(k)],['allersp_spca_common_groupcond_tworow','.pdf']);
%     exportgraphics(gcf,config.save_filepath_dropbox);
    config.save_filepath_M = fullfile(cluster_dir,'ERSP_Plots',['Cluster_',num2str(k)],['allersp_spca_common_groupcond_tworow','.pdf']);
    exportgraphics(gcf,config.save_filepath_M);
end


