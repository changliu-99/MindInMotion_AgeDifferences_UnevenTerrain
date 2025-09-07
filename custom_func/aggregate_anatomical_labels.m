function [atlas_name] = aggregate_anatomical_labels(cluster_update,k)

    if ispc
        FIELDTRIP_PATH = 'M:\liu.chang1\fieldtrip-master';
    else
        FIELDTRIP_PATH = '/blue/dferris/liu.chang1/fieldtrip-master';
    end
    
    atlas = ft_read_atlas([FIELDTRIP_PATH filesep 'template' filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii']);
    %- extract centroid locations
    centroids = [cluster_update(k).all_diplocs];
    
    dipfit_roi = cat(1,centroids);
    %- params
    atlas_name = cell(size(dipfit_roi,1),2);
    for i = 2:length(dipfit_roi)
        cfg              = [];
        %- (01/19/2023), JS: Perhaps have the function loop through each dipoel
        %in the cluster, determine a location of interest, and do it that way?
        %Doesn't seem like passing in a set of points helps...
    %     cfg.roi = STUDY.cluster(i).preclust.preclustdata;
        cfg.roi        = dipfit_roi(i,:);
        cfg.output     = 'multiple';
        cfg.atlas      = atlas;
        cfg.inputcoord = 'mni';
        %- (01/19/2023), JS: Changing cfg.sphere from 1 to 3 (1mm vs 3mm)
        cfg.sphere = 3;
        labels = ft_volumelookup(cfg, atlas);
        if ~isempty(labels)
            % disp(labels.count(labels.count ~= 0))
            [val, indx] = max(labels.count);
            if strcmp(labels.name(indx),'no_label_found')
                sub_indx = find(labels.count ~= 0 & labels.count < val);
                if isempty(sub_indx)
                    atlas_name{i,1} = ['CLs ',num2str(i)];
                    atlas_name{i,2} = labels.name{indx};
                    continue;
                end
                atlas_name{i,1} = ['CLs ',num2str(i)];
                tmp = sprintf('(N=%i) %s',labels.count(sub_indx(1)),labels.name{sub_indx(1)});
                for ic_i = sub_indx(2:end)
                    tmp = [tmp,' & ', sprintf('(N=%i) %s',labels.count(ic_i),labels.name{ic_i})];
                end
                atlas_name{i,2} = tmp;
            else
                atlas_name{i,1} = ['CLs ',num2str(i)];
                atlas_name{i,2} = labels.name{indx};
            end
        else
            warning('No label found for cluster %i',i);
        end
    end

    fprintf('Cluster \t\t Label\n');
    fprintf('________________________\n');
    for i = 2:size(dipfit_roi,1)
        label = cellstr(atlas_name{i,2});
        cl =  cellstr(atlas_name{i,1});
        fprintf('%s\t\t%s\n',cl{:},label{:})
    end
end