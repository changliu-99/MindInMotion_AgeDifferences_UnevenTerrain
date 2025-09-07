function [specdata_all_2group,specdata_all_3group] = organize_specdata(specdata_rest,specdata_gait,STUDY,cluster_num,cond_idx)
    group_list = {STUDY.datasetinfo(:).group};    

    % assign sets to each group
    STUDY.cluster(cluster_num).group = group_list(STUDY.cluster(cluster_num).sets);

    STUDY.cluster(cluster_num).group(strcmp(STUDY.cluster(cluster_num).group,'H1000s')) = {1};
    STUDY.cluster(cluster_num).group(strcmp(STUDY.cluster(cluster_num).group,'H2000s')) = {2};
    STUDY.cluster(cluster_num).group(strcmp(STUDY.cluster(cluster_num).group,'H3000s')) = {3};
    STUDY.cluster(cluster_num).group = cellfun(@(x) [[]; x], STUDY.cluster(cluster_num).group);
    
    specdata_all_2group = cell(5,2);
    specdata_all_2group{1,1} = specdata_rest{1};
    specdata_all_2group{1,2} = [specdata_rest{2} specdata_rest{3}];
    for i = cond_idx
        specdata_all_2group{i+1,1} = specdata_gait{i}(:,STUDY.cluster(cluster_num).group == 1);
        specdata_all_2group{i+1,2} = specdata_gait{i}(:,STUDY.cluster(cluster_num).group > 1);
    end
    
    specdata_all_3group = cell(5,3);
    specdata_all_3group{1,1} = specdata_rest{1};
    specdata_all_3group{1,2} = specdata_rest{2};
    specdata_all_3group{1,3} = specdata_rest{3};
    for i = cond_idx
        specdata_all_3group{i+1,1} = specdata_gait{i}(:,STUDY.cluster(cluster_num).group == 1);
        specdata_all_3group{i+1,2} = specdata_gait{i}(:,STUDY.cluster(cluster_num).group == 2);
        specdata_all_3group{i+1,3} = specdata_gait{i}(:,STUDY.cluster(cluster_num).group == 3);
    end
end

    
    
    
