%% Figure Plot original spec_data and the flattened data side by side 
function plot_original_and_flattened_spec(cfg,cls)
    color = cfg.color;
    spec_data_original = cfg.spec_data_original;
    fooof_group_results_org = cfg.fooof_group_results_org;
    fooof_diff_store = cfg.fooof_diff_store;
    fooof_freq = cfg.fooof_freq;
    fooof_apfit_store = cfg.fooof_apfit_store;

   
    PSD_YLIM_ORG = [-30 -10];
    PSD_YLIM_FOOOF = [-2 12];
    for g = 1
        switch g
            case 1
                color_dark = [color.terrain]; 
        %         color_dark(color_dark<0) = 0;
                color_light = [color.terrain_shade];
        %         color_dark =color.terrain(2:end,:);
        %         color_light =color.terrain_shade(2:end,:);
                xtick_label_g = {'rest','flat','low','med','high'};
            case 2
                color_dark = COLOR_MAPS_SPEED; %color.speed;
                color_light = COLOR_MAPS_SPEED+0.15; %color.speed_shade;
                xtick_label_g = {'0.25','0.5','0.75','1.0'};
        end
        % figure('color','white','position',[100 300 600 1000],'Renderer','Painters');
        j = 1;
        for k = cls
            if length(spec_data_original{g}{k})==4
                color_dark = color_dark(2:end,:);
                color_light = color_light(2:end,:);
                xtick_label_g = xtick_label_g{2:end};
            end
            brain_area = cfg.brain_area;
                switch brain_area
                    case 'SMA'
                        ylim_intv_orig = [-30 -10];
                        ylim_intv_fooof = [-1 12];
                    case 'PP'
                        ylim_intv_orig = [-30 -5];
                        ylim_intv_fooof = [-1 12];
                    case 'Precuneus'
                        ylim_intv_orig = [-30 -10];
                        ylim_intv_fooof = [-1 11];
                    case 'Occipital'
                        ylim_intv_orig = [-30 -8];
                        ylim_intv_fooof = [-1 11];
                    case 'SuppMotor'
                        ylim_intv_orig = [-30 -12];
                        ylim_intv_fooof = [-1 5];
                    case {'','Cingulate'} 
                        ylim_intv_orig = [-30 -10];
                        ylim_intv_fooof = [-1 5];
                end
            % hardcode to preset the axis limit
            %## set color limits
            data_1 = spec_data_original{g}{k}{1}';
            data_2 = cat(2,fooof_diff_store{g}{k}{:});
                       
%             ylim_intv_orig = [min(mean(data_1))*1.1 max(mean(data_1))*0.7];%PSD_YLIM_ORG; %[-30,-10]; %zeros(1,2);
%             ylim_intv_fooof = [-1 max(mean(data_2))]*2; %[-2,6]; %zeros(1,2);
            if mod(j,2) == 1 
                figure('color','white','position',[100 300 1020 180],'Renderer','Painters');
                j = 1;
            end

            % -------------- original PSD -----------
            subplot(1,4,j)
            title('Young')
            for i = 1:length(spec_data_original{g}{k})
                data = spec_data_original{g}{k}{i,1}';
                JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));% input need to be a row vector
            end
            for i = 1:length(spec_data_original{g}{k})
                data = spec_data_original{g}{k}{i,1}';
                plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
            end  
            % plot the aperiodic line
            for i = 1:length(spec_data_original{g}{k})
                aperiodic_fit = fooof_apfit_store{g}{k}{i}';
                plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',0.5);
            end
            ylim(ylim_intv_orig);
            ax = gca;       
            axsignif = highlight_CL(ax, fooof_freq, [], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('');
            ylabel('10*log(PSD)');
            set(gca,'fontsize',12);
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
        %     title(['Cluster ',num2str(k)]);
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            subplot(1,4,j+1)
            title('Older')
            for i = 1:length(spec_data_original{g}{k})
                data = spec_data_original{g}{k}{i,2}';
                JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));% input need to be a row vector
            end
            for i = 1:length(spec_data_original{g}{k})
                data = spec_data_original{g}{k}{i,2}';
                plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
            end  
            % plot the aperiodic line
            for i = 1:length(spec_data_original{g}{k})
                aperiodic_fit = fooof_apfit_store{g}{k}{i,2}';
                plot(fooof_freq,mean(aperiodic_fit),'color',color_dark(i,:),'linestyle','--','linewidth',0.5);
            end
            ylim(ylim_intv_orig);
            ax = gca;       
            axsignif = highlight_CL(ax, fooof_freq, [], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('');
            ylabel('');
            set(gca,'fontsize',12);
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
        %     title(['Cluster ',num2str(k)]);
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            % -------------- Flattened PSD -----------
            subplot(1,4,j+2)
            title('Young');
            for i = 1:length(fooof_diff_store{g}{k})
                data = fooof_diff_store{g}{k}{i,1}';
                JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));% input need to be a row vector
            end
            for i = 1:length(fooof_diff_store{g}{k})
                data = fooof_diff_store{g}{k}{i,1}';
                plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
            end  
            ylim(ylim_intv_fooof);
            ax = gca;       
            axsignif = highlight_CL(ax, fooof_freq,[], 'background', '');
            xlim([4 40]);
            plot([0 40],[0 0],'-','color',[0.5 0.5 0.5]);
            ylabel('');
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            subplot(1,4,j+3)
            title('Old');
            for i = 1:length(fooof_diff_store{g}{k})
                data = fooof_diff_store{g}{k}{i,2}';
                JackKnife_sung(fooof_freq,mean(data),[mean(data)-std(data)/sqrt(size(data,1))],[mean(data)+std(data)/sqrt(size(data,1))],...
                    color_dark(i,:),color_light(i,:));% input need to be a row vector
            end
            for i = 1:length(fooof_diff_store{g}{k})
                data = fooof_diff_store{g}{k}{i,2}';
                plot(fooof_freq,mean(data),'color',color_dark(i,:),'linewidth',1);
            end  
            ylim(ylim_intv_fooof);
            ax = gca;       
            axsignif = highlight_CL(ax, fooof_freq,[], 'background', '');
            xlim([4 40]);
            plot([0 40],[0 0],'-','color',[0.5 0.5 0.5]);
            ylabel('');
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            j = j + 4;

%             if g == 1
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'.pdf']),'ContentType','vector','Resolution',300)
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'.jpg']),'ContentType','vector','Resolution',300)
%             elseif g == 2
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'_speed.pdf']),'ContentType','vector','Resolution',300)
%             end
        end
    end
end

