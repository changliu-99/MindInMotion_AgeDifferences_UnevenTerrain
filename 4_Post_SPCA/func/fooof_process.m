function [fooof_outcome] = fooof_process(specdata_all_2group,TMP_STUDY,cluster_update,k,g,cfg)
    
    f_range = cfg.f_range;
    settings = cfg.settings;
    SUB_GROUP_FNAME_REGEX = [];
    specfreqs_rest = cfg.specfreqs_rest;
    
    theta_band = [4, 8];
    alpha_band = [8 12];
    beta_band  = [12 30];

    %% Run fooof
    % Input should be in linear spacing   
    i_ind = 0;
    %- get subjects in cluster
    s_chars = {TMP_STUDY.datasetinfo(cluster_update(k).sets).subject};
    
    cl_chars = s_chars;
    cl_inds = cluster_update(k).sets;
    cl_comps = cluster_update(k).comps;
    cl_speeds = zeros(length(cl_chars),1);
    
    num_young = size(specdata_all_2group{1,1},2);
%         for i = 1:length(cl_speeds)
%             ind = cellfun(@(x) x == categorical(cl_chars(i)),speed_alleeg(:,1));
%             cl_speeds(i) = speed_alleeg{ind,2};
%         end
    for group = 1:size(specdata_all_2group,2) % in case there is young and old adult group
        for cond = 1:size(specdata_all_2group,1) % different level of terrains
            specdata_nolog = 10.^(specdata_all_2group{cond,group}/10);
            % Run FOOOF
            return_model = true;
            for i = 1:size(specdata_all_2group{cond,group},2)
                if ~sum(isnan(specdata_nolog(:,i)))
                    fooof_results{g}{cond,group}{i} = fooof(specfreqs_rest, specdata_nolog(:,i), f_range, settings, return_model);
    %                     fooof_group_results_org{g}{k}(i_ind + i).subjects = s_chars(g_inds);
                    fooof_group_results_org{g}{k}(i_ind + i).speed_ms = cl_speeds(i+(group-1)*num_young);
                    fooof_group_results_org{g}{k}(i_ind + i).subID = cl_inds(i+(group-1)*num_young); %cluster_update(k).sets(i);
                    fooof_group_results_org{g}{k}(i_ind + i).sub_char = cl_chars{i+(group-1)*num_young}; %cluster_update(k).sets(i);
                    fooof_group_results_org{g}{k}(i_ind + i).compID = cl_comps(i+(group-1)*num_young); %cluster_update(k).comps(i);
                    %-
                    fooof_group_results_org{g}{k}(i_ind + i).study = g;%1 = terrain, 2 = speed
                    fooof_group_results_org{g}{k}(i_ind + i).cond = cond;
                    fooof_group_results_org{g}{k}(i_ind + i).group = group;
                    fooof_group_results_org{g}{k}(i_ind + i).cluster = k;
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_exp = fooof_results{g}{cond,group}{i}.aperiodic_params(2);
                    fooof_group_results_org{g}{k}(i_ind + i).aperiodic_offset = fooof_results{g}{cond,group}{i}.aperiodic_params(1);
                    fooof_group_results_org{g}{k}(i_ind + i).central_freq = fooof_results{g}{cond,group}{i}.peak_params(:,1);
                    fooof_group_results_org{g}{k}(i_ind + i).power = fooof_results{g}{cond,group}{i}.peak_params(:,2);
                    fooof_group_results_org{g}{k}(i_ind + i).r_squared = fooof_results{g}{cond,group}{i}.r_squared;
                    %------------ Compute average power after flatten curve
    %                     fooof_diff = fooof_results{g}{cond,group}{i}.power_spectrum - fooof_results{g}{cond,group}{i}.ap_fit;
                    % Super important, the output is already logged, the
                    % only difference is the magnitude by 10
                    fooof_diff = 10*(fooof_results{g}{cond,group}{i}.power_spectrum) - 10*(fooof_results{g}{cond,group}{i}.ap_fit);
                    fooof_freq = fooof_results{g}{cond,group}{i}.freqs;
                    fooof_group_results_org{g}{k}(i_ind + i).theta_avg_power = mean(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_avg_power = mean(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).beta_avg_power = mean(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));

                    % find the min and max power
                    fooof_group_results_org{g}{k}(i_ind + i).theta_min_power = min(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_min_power = min(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).beta_min_power = min(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));
                    
                    fooof_group_results_org{g}{k}(i_ind + i).theta_max_power = max(fooof_diff(fooof_freq >= theta_band(1) & fooof_freq < theta_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_max_power = max(fooof_diff(fooof_freq >= alpha_band(1) & fooof_freq < alpha_band(2)));
                    fooof_group_results_org{g}{k}(i_ind + i).beta_max_power = max(fooof_diff(fooof_freq >= beta_band(1) & fooof_freq < beta_band(2)));                    
                    
                    % find at which frequency the min and max power occur
                    fooof_group_results_org{g}{k}(i_ind + i).theta_min_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).theta_min_power);
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_min_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).alpha_min_power);
                    fooof_group_results_org{g}{k}(i_ind + i).beta_min_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).beta_min_power);
                    fooof_group_results_org{g}{k}(i_ind + i).theta_max_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).theta_max_power);
                    fooof_group_results_org{g}{k}(i_ind + i).alpha_max_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).alpha_max_power);
                    fooof_group_results_org{g}{k}(i_ind + i).beta_max_power_freq = fooof_freq(fooof_diff == fooof_group_results_org{g}{k}(i_ind + i).beta_max_power);
                    
                    % data structure needs to be freq x subject
                    fooof_diff_store{g}{k}{cond,group}(:,i) = fooof_diff';
                    fooof_apfit_store{g}{k}{cond,group}(:,i) = 10*(fooof_results{g}{cond,group}{i}.ap_fit);

                    % - store original spec data
                    spec_data_original{g}{k}{cond,group} = specdata_all_2group{cond,group}(specfreqs_rest >= f_range(1) & specfreqs_rest <= f_range(2),:);
                end
            end
            i_ind = i_ind + size(specdata_all_2group{cond,group},2);
        end
    end
    fooof_outcome.fooof_group_results_org = fooof_group_results_org;
    fooof_outcome.fooof_diff_store = fooof_diff_store;
    fooof_outcome.fooof_apfit_store = fooof_apfit_store;
    fooof_outcome.spec_data_original = spec_data_original;
end