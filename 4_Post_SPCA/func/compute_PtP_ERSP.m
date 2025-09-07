% Compute the Peak to Peak ERSP fluctuation
% Chang Liu
% 20240208
function p2p_ERSP_table_all = compute_PtP_ERSP(allersp,alltimes,allfreqs,config)
    % input: allersp = freqs x time x subject;
    % setup the output struct
    p2p_ERSP_table_all =table;
    
    alpha = 0.05;
    warpingvalues = config.warpingvalues;
    baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));  
    YlimMax = config.YlimMax;
        
    % group ERSP into theta,alpha, beta band
    theta_band = allfreqs <= 8 & allfreqs > 3;
    alpha_band = allfreqs <= 13 & allfreqs > 8;
    beta_band = allfreqs <= 30 & allfreqs > 13;
    gamma_band = allfreqs <= 50 & allfreqs > 30;
    
    var_name = {theta_band;alpha_band;beta_band;gamma_band};
    title_keyword = {'theta','alpha','beta','gamma'};
    % plot
%     fig = figure('color','white','position',[200 200 800 200]);
    for band = 1:4
%         subplot(1,4,band)
        for j = 1:4
            % organize to table
            temp = squeeze(mean(allersp{j}(var_name{band},:,:),1));
            p2p_ERSP = squeeze(peak2peak(temp(baseidx,:)));
            p2p_ERSP_table_tmp = table(repmat(title_keyword(band),length(p2p_ERSP),1),repmat(j+1,length(p2p_ERSP),1),p2p_ERSP',...
                'VariableNames',{'band','cond','p2p_ersp'});
            p2p_ERSP_table_all = vertcat(p2p_ERSP_table_all, p2p_ERSP_table_tmp);
            allersp_to_plot = median(mean(allersp{j}(var_name{band},:,:),1),3);%take mean across freqs and then take mean across participants
%             plot(alltimes(:,baseidx),allersp_to_plot(1,baseidx,:),'color',config.color.terrain(j+1,:),'linewidth',2);
%             hold on;
            clear temp
        end
%         plot(alltimes(:,baseidx),zeros(1,length(baseidx)),'--','color','black');
%         set(gca,'xlim',[warpingvalues(1) warpingvalues(end)])
%         
%         if band == 1
%             ylabel(sprintf('dB'),'fontsize',10,'fontweight','bold');
%         else
%             ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
%         end
%         xlabel('','Fontsize',12);
% %             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
%         title(title_keyword{band});
%         set(gca,'xtick',warpingvalues,'xticklabel',{'RHS','LTO','LHS','RTO','RHS'});
%         xtickangle(45)
%         ax = gca; ax.XAxis.FontSize = 8;
        
    end
    if config.save_file
        saveas(fig,config.save_filepath);
        close all
    end
end