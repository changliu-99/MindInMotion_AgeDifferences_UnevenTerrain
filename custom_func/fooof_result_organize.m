function [psd_feature_table,psd_feature_table_unstack,psd_feature_table_simplify] = fooof_result_organize(fooof_group_results_org,DESIGN_I,save_filename)
    % function that organize the fooof table together
    
    theta_band = [4, 8];
    alpha_band = [8 12];
    beta_band  = [12 30];
    cond_terrains = {'flat','low','med','high','rest'};
    cond_speeds = {'0p25','0p5','0p75','1p0','rest'};

    %% Determine peak occurs at different band, don't think I am using this
    % Looks like 
    % Define frequency bands of interest
    for g = 1
        for k = 3:length(fooof_group_results_org{g})
            for i = 1:length(fooof_group_results_org{g}{k})
                if ~isempty(fooof_group_results_org{g}{k}(i).central_freq)
                    fooof_group_results_org{g}{k}(i).theta = [];
                    fooof_group_results_org{g}{k}(i).alpha = [];
                    fooof_group_results_org{g}{k}(i).beta = [];
                    for j = 1:length(fooof_group_results_org{g}{k}(i).central_freq)
                        cf = fooof_group_results_org{g}{k}(i).central_freq(j);
                        if cf > theta_band(1) & cf <= theta_band(2)
                            fooof_group_results_org{g}{k}(i).theta = [fooof_group_results_org{g}{k}(i).theta; cf fooof_group_results_org{g}{k}(i).power(j)];
                        elseif cf > alpha_band(1) & cf <= alpha_band(2)
                            fooof_group_results_org{g}{k}(i).alpha = [fooof_group_results_org{g}{k}(i).alpha; cf fooof_group_results_org{g}{k}(i).power(j)];
                        elseif cf > beta_band(1) & cf <= beta_band(2)
                            fooof_group_results_org{g}{k}(i).beta = [fooof_group_results_org{g}{k}(i).beta; cf fooof_group_results_org{g}{k}(i).power(j)];
                        end
                    end
                    if length(fooof_group_results_org{g}{k}(i).theta) > 1
                        [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).theta(:,1)-6));
                        temp_power = fooof_group_results_org{g}{k}(i).theta(:,2);
                        fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).theta(indx,1) temp_power(indx)];
                    end
                    if length(fooof_group_results_org{g}{k}(i).alpha) > 1
                        [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).alpha(:,1)-10));
                        temp_power = fooof_group_results_org{g}{k}(i).alpha(:,2);
                        fooof_group_results_org{g}{k}(i).alpha_center = [fooof_group_results_org{g}{k}(i).alpha(indx,1) temp_power(indx)];
                    end
                    if length(fooof_group_results_org{g}{k}(i).beta) > 1
                        [~,indx] = min(abs(fooof_group_results_org{g}{k}(i).beta(:,1)-20));
                        temp_power = fooof_group_results_org{g}{k}(i).beta(:,2);
                        fooof_group_results_org{g}{k}(i).beta_center = [fooof_group_results_org{g}{k}(i).beta(indx,1) temp_power(indx)];
                    end
                end
            end
        end
    end
    %## SAVE DATA
    % save([save_dir filesep 'fooof_results_summary.mat'],'fooof_group_results_org');
    % save([save_dir filesep 'fooof_diff_store.mat'],'fooof_diff_store');
    % save([save_dir filesep 'fooof_apfit_store.mat'],'fooof_apfit_store');
    % save([save_dir filesep 'spec_data_original.mat'],'spec_data_original');
    % 
    %% Create table from group results, take mean across participants ICs    
    psd_feature_table = table;
    for g = DESIGN_I
        for k = 3:length(fooof_group_results_org{g})
            disp(k)
            temp_table_C1 = fooof_group_results_org{g}{k};
            for i = 1:length(temp_table_C1)
                temp_table_C1(i).alpha = [];
                temp_table_C1(i).beta = [];
                temp_table_C1(i).alpha_center = [];
                temp_table_C1(i).beta_center = [];
                if isempty(temp_table_C1(i).alpha)
                    temp_table_C1(i).alpha = [NaN, NaN];
                end
                if isempty(temp_table_C1(i).beta)
                    temp_table_C1(i).beta = [NaN, NaN];
                end
                if isempty(temp_table_C1(i).alpha_center)
                    temp_table_C1(i).alpha_center = [NaN, NaN];
                end

                if isempty(temp_table_C1(i).beta_center)
                    temp_table_C1(i).beta_center = [NaN, NaN];
                end

                [~,idx_a] = max(temp_table_C1(i).alpha(:,2));
                [~,idx_b] = max(temp_table_C1(i).beta(:,2));
                
                try
                    psd_feature_table = vertcat(psd_feature_table,table(temp_table_C1(i).subID,temp_table_C1(i).compID,temp_table_C1(i).study,...
                        temp_table_C1(i).cond,temp_table_C1(i).group,temp_table_C1(i).cluster,cellstr(temp_table_C1(i).sub_char),temp_table_C1(i).speed_ms,temp_table_C1(i).aperiodic_exp,temp_table_C1(i).aperiodic_offset,...
                        temp_table_C1(i).r_squared,temp_table_C1(i).alpha(idx_a,1),temp_table_C1(i).alpha(idx_a,2),temp_table_C1(i).beta(idx_b,1),temp_table_C1(i).beta(idx_b,2),...
                        temp_table_C1(i).alpha_center(1,1),temp_table_C1(i).alpha_center(1,2),temp_table_C1(i).beta_center(1,1),temp_table_C1(i).beta_center(1,2),...
                        temp_table_C1(i).theta_avg_power(1,1),temp_table_C1(i).alpha_avg_power(1,1),temp_table_C1(i).beta_avg_power(1,1),...
                        temp_table_C1(i).theta_min_power(1,1),temp_table_C1(i).alpha_min_power(1,1),temp_table_C1(i).beta_min_power(1,1),...
                        temp_table_C1(i).theta_max_power(1,1),temp_table_C1(i).alpha_max_power(1,1),temp_table_C1(i).beta_max_power(1,1),...
                        temp_table_C1(i).theta_min_power_freq(1,1),temp_table_C1(i).alpha_min_power_freq(1,1),temp_table_C1(i).beta_min_power_freq(1,1),...
                        temp_table_C1(i).theta_max_power_freq(1,1),temp_table_C1(i).alpha_max_power_freq(1,1),temp_table_C1(i).beta_max_power_freq(1,1),...
                        'VariableNames',...
                        {'subID','compID','study','cond','group','cluster','subj_char','speed_ms','aperiodic_exp','aperiodic_offset','r_squared','alpha_cf','alpha_p',...
                        'beta_cf','beta_p','alpha_center','alpha_centerP','beta_center','beta_centerP',...
                        'theta_avg_power','alpha_avg_power','beta_avg_power',...
                        'theta_min_power','alpha_min_power','beta_min_power','theta_max_power','alpha_max_power','beta_max_power',...
                        'theta_min_power_freq','alpha_min_power_freq','beta_min_power_freq','theta_max_power_freq','alpha_max_power_freq','beta_max_power_freq'}));
                catch
                    disp('skipping one row');
                end
            end
        end
    end
    psd_feature_table.subID = categorical(psd_feature_table.subID);
    psd_feature_table.cond = categorical(psd_feature_table.cond);
    psd_feature_table.cluster = categorical(psd_feature_table.cluster);
    psd_feature_table.study = categorical(psd_feature_table.study);

    % C1_table.study = categorical(C1_table.study);

%     grp_C1_table = grpstats(psd_feature_table,["subID","study","cond","cluster"],'nanmedian','DataVars',["r_squared","alpha_cf","alpha_p",...
%                 "beta_cf","beta_p","alpha_center","alpha_centerP","beta_center","beta_centerP","theta_avg_power","alpha_avg_power","beta_avg_power"]);
%     grp_C1_table.Properties.VariableNames = {'subID','study','cond','cluster','GroupCount','med_r_squared','med_alpha_cf','med_alpha_p',...
%                 'med_beta_cf','med_beta_p','med_alpha_center','med_alpha_centerP','med_beta_center','med_beta_centerP',...
%                 'med_theta_avg_power','med_alpha_avg_power','med_beta_avg_power'};

    % -- Reorganize table
    psd_feature_table_simplify = table(psd_feature_table.subID,psd_feature_table.subj_char,psd_feature_table.study,psd_feature_table.group,...
        psd_feature_table.cond,psd_feature_table.cluster,psd_feature_table.theta_avg_power,...
        psd_feature_table.alpha_avg_power,psd_feature_table.beta_avg_power);
    psd_feature_table_simplify.Properties.VariableNames = {'subID','subjectName','study','group','cond','cluster','theta_avg_power','alpha_avg_power','beta_avg_power'};
    psd_feature_table_unstack = unstack(psd_feature_table_simplify,{'theta_avg_power','alpha_avg_power','beta_avg_power'},'cond','VariableNamingRule','preserve');

    % -- SAVE table
    writetable(psd_feature_table,[save_filename.save_dir filesep [save_filename.xlsx_name,'.xlsx']]);
    save([save_filename.save_dir filesep [save_filename.psd_table_name,'.mat']],'psd_feature_table');
    save([save_filename.save_dir filesep [save_filename.psd_table_unstack_name,'.mat']],'psd_feature_table_unstack');
    
    % -- SAVE simplified table that only contained relevant variables - update 20240425
    writetable(psd_feature_table_simplify,[save_filename.save_dir filesep [save_filename.xlsx_name,'_simplify.xlsx']]);
    save([save_filename.save_dir filesep [save_filename.psd_table_name,'_simplify.mat']],'psd_feature_table_simplify');
    
