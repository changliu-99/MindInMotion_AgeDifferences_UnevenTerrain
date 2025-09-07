function [comps_out,main_cl_inds,outlier_cl_inds,valid_clusters,main_cl_anat,nonzero_cl_inds] = eeglab_get_cluster_comps_CL(STUDY,varargin)
%EEGLAB_GET_CLUSTER_COMPS Summary of this function goes here
%   Function generates a matrix (NxM, N = number of clusters, M = number of
%   subjects) of DOUBLES. Each element in the matrix is the component
%   number for each i'th subject in STUDY.datasetinfo/ALLEEG.
%   IN: 
%   OUT: 
%   IMPORTANT: 
% Code Designer: Jacob Salminen (11/25/2022)
% Chang Liu: adapt for both young and old adult group
%## TIME
tic
%## (PARSER) DEFINE DEFAULTS

%## (PARSER) Define Parser
p = inputParser;
%## (PARSER) REQUIRED
addRequired(p,'STUDY',@isstruct)
%## (PARSER) OPTIONAL
%## (PARSER) PARAMETER
parse(p,STUDY,varargin{:});
%## (PARSER) SET DEFAULTS
%- OPTIONALS
%- PARAMETER
%## (PARSER) MAKE DIRS
%% ===================================================================== %%
main_cl_anat = cell(length(length(STUDY.cluster)),1);
%## Extract components for each cluster & subject
%- PARAMS
subj_inds = unique(STUDY.cluster(1).sets);
comps_out = zeros(length(STUDY.cluster),length(subj_inds));
main_cl_inds = zeros(1,length(STUDY.cluster));
outlier_cl_inds = zeros(1,length(STUDY.cluster));
%- loop through all clusters but the parent cluster
main_cl_inds(1) = 1;
for clus_i = 2:length(STUDY.cluster)
    %- if cluster name contains the word 'Outlier' or 'Reject' then the cluster is an
    %outlier
    chk = isempty(regexpi(STUDY.cluster(clus_i).name,'outlier')) && isempty(regexpi(STUDY.cluster(clus_i).name,'reject'));
    if chk
        sets_i = STUDY.cluster(clus_i).sets;
        main_cl_inds(clus_i) = 1;
        for j = 1:length(sets_i)
            comps_out(clus_i,sets_i(j)) = STUDY.cluster(clus_i).comps(j);
        end
    else
        outlier_cl_inds(clus_i) = 1;
    end
end
tmp = cellfun(@length,cellfun(@unique,{STUDY.cluster.sets},'UniformOutput',false),'UniformOutput',false);
tmp = cell2mat(tmp);
%- clusters that have comps
nonzero_cl_inds = zeros(1,length(STUDY.cluster));
for clus_i = 2:length(STUDY.cluster)
    nonzero_cl_inds(clus_i) = ~isempty(STUDY.cluster(clus_i).comps);
end
nonzero_cl_inds = find(nonzero_cl_inds);
%- clusters with >50% of the subjects

% -- added by Chang
num_young = sum(contains(STUDY.subject,'H1'));
num_old = sum(contains(STUDY.subject,'H2')) + sum(contains(STUDY.subject,'H3'));
group_list = {STUDY.datasetinfo(:).group};

j = 1;
for i = 3:length(STUDY.cluster)
    STUDY.cluster(i).group = group_list(STUDY.cluster(i).sets);
    STUDY.cluster(i).group(strcmp(STUDY.cluster(i).group,'H1000s')) = {1};
    STUDY.cluster(i).group(strcmp(STUDY.cluster(i).group,'H2000s')) = {2};
    STUDY.cluster(i).group(strcmp(STUDY.cluster(i).group,'H3000s')) = {2};
    STUDY.cluster(i).group = cellfun(@(x) [[]; x], STUDY.cluster(i).group);
    group_store = STUDY.cluster(i).group;
    
    if sum(group_store==1) >=0.5*((num_young)) & sum(group_store==2) >=0.5*((num_old)) & i < sum(main_cl_inds)+2
        valid_clusters(j) = i;
        j = j + 1;
    end
    
end
%- remove all outlier cluster components
% comps_out = comps_out(logical(main_cl_inds),:);
main_cl_inds = find(main_cl_inds);
outlier_cl_inds = find(outlier_cl_inds);
if isfield(STUDY.cluster,'analabel')
    main_cl_anat = {STUDY.cluster.analabel};
end
if all(comps_out==0)
    error('STUDY cluster information not generated');
end
end

