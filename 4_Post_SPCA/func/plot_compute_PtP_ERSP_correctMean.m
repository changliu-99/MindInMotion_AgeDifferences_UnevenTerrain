% Compute the Peak to Peak ERSP fluctuation with mean correction for ERSP
% Chang Liu
% 20240208
function p2p_ERSP_table_all = plot_compute_PtP_ERSP_correctMean(allersp,alltimes,allfreqs,config)
    % input: allersp = freqs x time x subject;
    % setup the output struct
    p2p_ERSP_table_all =table;
    
    alpha = 0.05;
    warpingvalues = config.warpingvalues;
    baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));  
    YlimMax = config.YlimMax;    
    color_dark = [config.color.terrain]; 
    color_light = [config.color.terrain_shade];
    
    % correct mean by subtracting the mean from each participant
    for g = 1:size(allersp,2)
        for j = 1: length(allersp)
            erspdata = allersp{j,g}(:,:,:);
            baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));                
            baseline_allcomp = mean(erspdata(:,baseidx,:),2); % mean power for each person
            baseline = mean(baseline_allcomp,3);%mean power across participant
            allersp_sub_mean{j,g} = allersp{j,g}(:,:,:)-repmat(baseline_allcomp,1,length(alltimes));% subtract baseline for each person
    %             cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline,1,length(alltimes)),3);
        end
    end
    % group ERSP into theta,alpha, beta band
    theta_band = allfreqs <= 8 & allfreqs > 3;
    alpha_band = allfreqs <= 13 & allfreqs > 8;
    beta_band = allfreqs <= 30 & allfreqs > 13;
    gamma_band = allfreqs <= 50 & allfreqs > 30;
    
    var_name = {theta_band;alpha_band;beta_band;gamma_band};
    title_keyword = {'\theta','\alpha','\beta','gamma'};
    % plot
    %{
    fig = figure('color','white','position',[200 200 300 300]);hold on;
    p = 1;
    for group = 1:2
        for band = 1:3
            subplot(2,3,p)
            for j = 1:4
                % organize to table
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = median(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                std_to_plot = std(data(:,baseidx,:),0,3);
%                 plot(alltimes(:,baseidx),allersp_to_plot(1,baseidx,:),'color',config.color.terrain(j+1,:),'linewidth',2);
                JackKnife_sung(alltimes(:,baseidx),allersp_to_plot,...
                    [allersp_to_plot-std_to_plot/sqrt(size(data,3))],[allersp_to_plot+std_to_plot/sqrt(size(data,3))],...
                    color_dark(j+1,:),color_light(j+1,:));
                hold on;
            end
            for j = 1:4
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = median(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),allersp_to_plot,'color',color_dark(j+1,:),'linewidth',1);
            end
            plot(alltimes(:,baseidx),zeros(1,length(baseidx)),'-','color','black');
            set(gca,'xlim',[warpingvalues(1) warpingvalues(end)])

            if band == 1
                ylabel(sprintf('Power (dB)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
    %             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
            if group == 1;title(title_keyword{band});else;title('');end
            set(gca,'xtick',warpingvalues,'xticklabel',{'','','','',''});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;ax.YAxis.FontSize = 10;
        
            ylim([-config.climMat_max config.climMat_max]);
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4) alltimes(baseidx(end))],{'k--' ,'k--', 'k--', 'k--','k--'});

            p = p + 1;
        end
    end
    %}
    fig = figure('color','white','position',[200 200 300 600]);hold on;
    p = 1;
    for band = 1:3
        for group = 1:2
            subplot(3,2,p)
            for j = 1:4
                % organize to table
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = mean(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                std_to_plot = std(data(:,baseidx,:),0,3);
%                 plot(alltimes(:,baseidx),allersp_to_plot(1,baseidx,:),'color',config.color.terrain(j+1,:),'linewidth',2);
                JackKnife_sung(alltimes(:,baseidx),allersp_to_plot,...
                    [allersp_to_plot-std_to_plot/sqrt(size(data,3))],[allersp_to_plot+std_to_plot/sqrt(size(data,3))],...
                    color_dark(j+1,:),color_light(j+1,:));
                hold on;
            end
            for j = 1:4
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = mean(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),allersp_to_plot,'color',color_dark(j+1,:),'linewidth',1);
            end
            plot(alltimes(:,baseidx),zeros(1,length(baseidx)),'-','color','black');
            set(gca,'xlim',[warpingvalues(1) warpingvalues(end)])

            if group == 1
                ylabel(sprintf('%s Power (dB)',title_keyword{band}),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
    %             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
%             title(title_keyword{band});else;title('');end
            
            set(gca,'xtick',warpingvalues,'xticklabel',{'','','','',''});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;ax.YAxis.FontSize = 10;
        
            ylim([-config.climMat_max config.climMat_max]);
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4) alltimes(baseidx(end))],{'k--' ,'k--', 'k--', 'k--','k--'});

            p = p + 1;
        end
    end
    
    if config.save_file
        exportgraphics(fig,config.save_filepath);
        exportgraphics(fig,config.save_pdf_filepath);
        close all
    end
    %% Plot individual trajectory
    %{
    p = 1;
    fig = figure('color','white','position',[200 200 500 300]);hold on;
    for group = 1:2
        for band = 1:3
            subplot(2,3,p);hold on;
            for j = config.cond
                % organize to table
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = data(:,baseidx,:);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),squeeze(allersp_to_plot),'color',color_light(j+1,:),'linewidth',2);
                
                hold on;
            end
            for j = config.cond
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = mean(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),allersp_to_plot,'color',color_dark(j+1,:),'linewidth',1);
            end
            plot(alltimes(:,baseidx),zeros(1,length(baseidx)),'-','color','black');
            set(gca,'xlim',[warpingvalues(1) warpingvalues(end)])

            if band == 1
                ylabel(sprintf('Power (dB)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
    %             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
            title(title_keyword{band});
            set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;ax.YAxis.FontSize = 10;
        
            ylim([-3 3]);
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});

            p = p + 1;
        end
    end
    %}
    p = 1;
    fig = figure('color','white','position',[200 200 300 500]);hold on;
    for band = 1:3
        for group = 1:2
            subplot(3,2,p);hold on;
            for j = config.cond
                % organize to table
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = data(:,baseidx,:);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),squeeze(allersp_to_plot),'color',color_light(j+1,:),'linewidth',2);
                
                hold on;
            end
            for j = config.cond
                temp = squeeze(mean(allersp_sub_mean{j,group}(var_name{band},:,:),1));
                data = mean(allersp_sub_mean{j,group}(var_name{band},:,:),1);
                allersp_to_plot = mean(data(:,baseidx,:),3);%take mean across freqs and then take mean across participants
                plot(alltimes(:,baseidx),allersp_to_plot,'color',color_dark(j+1,:),'linewidth',1);
            end
            plot(alltimes(:,baseidx),zeros(1,length(baseidx)),'-','color','black');
            set(gca,'xlim',[warpingvalues(1) warpingvalues(end)])

            if band == 1
                ylabel(sprintf('Power (dB)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
    %             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
            title(title_keyword{band});
            set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;ax.YAxis.FontSize = 10;
        
            ylim([-3 3]);
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});

            p = p + 1;
        end
    end
    %% % ORGANIZE into a table
    for j = 1:4
        allersp_all_group_submean{j,1} = cat(3,allersp_sub_mean{j,1},allersp_sub_mean{j,2});
    end
    for band = 1:4
%         subplot(1,4,band)
        for j = 1:4
            % organize to table
            temp = squeeze(mean(allersp_all_group_submean{j}(var_name{band},:,:),1));
            p2p_ERSP = squeeze(peak2peak(temp(baseidx,:)));
            min_ERSP = squeeze(min(temp(baseidx,:)));
            max_ERSP = squeeze(max(temp(baseidx,:)));
            
            p2p_ERSP_table_tmp = table(repmat(title_keyword(band),length(p2p_ERSP),1),repmat(j+1,length(p2p_ERSP),1),p2p_ERSP',min_ERSP',max_ERSP',...
                'VariableNames',{'band','cond','p2p_ersp','min_ersp','max_ersp'});
            p2p_ERSP_table_all = vertcat(p2p_ERSP_table_all, p2p_ERSP_table_tmp);
%             plot(alltimes(:,baseidx),allersp_to_plot(1,baseidx,:),'color',config.color.terrain(j+1,:),'linewidth',2);
%             hold on;
            clear temp
        end
    end
    
    
end