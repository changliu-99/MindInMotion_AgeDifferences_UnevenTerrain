%% Figure Plot original spec_data and the flattened data side by side 
function plot_flattened_spec_stats(cfg,cls)
    color = cfg.color;
    spec_data_original = cfg.spec_data_original;
    fooof_group_results_org = cfg.fooof_group_results_org;
    fooof_diff_store = cfg.fooof_diff_store;
    fooof_freq = cfg.fooof_freq;
    fooof_apfit_store = cfg.fooof_apfit_store;

    pcond = cfg.stats.pcond;
    pgroup = cfg.stats.pgroup;
    pinter = cfg.stats.pinter;
    
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
%                         ylim_intv_fooof = [-1 12];
                        ylim_intv_fooof = [-1 6];
                    case 'PP'
                        ylim_intv_orig = [-30 -5];
%                         ylim_intv_fooof = [-1 12];
                        ylim_intv_fooof = [-1 8];
                    case 'Precuneus'
                        ylim_intv_orig = [-30 -10];
%                         ylim_intv_fooof = [-1 11];
                        ylim_intv_fooof = [-1 8];
                    case 'Occipital'
                        ylim_intv_orig = [-30 -8];
%                         ylim_intv_fooof = [-1 11];
                        ylim_intv_fooof = [-1 8];
                    case 'SuppMotor'
                        ylim_intv_orig = [-30 -12];
%                         ylim_intv_fooof = [-1 5];
                        ylim_intv_fooof = [-1 4];
                    case {'','Cingulate'} 
                        ylim_intv_orig = [-30 -10];
%                         ylim_intv_fooof = [-1 5];
                        ylim_intv_fooof = [-1 5];
                end
            % hardcode to preset the axis limit
            %## set color limits
            data_1 = spec_data_original{g}{k}{1}';
            data_2 = cat(2,fooof_diff_store{g}{k}{:});
                       
%             ylim_intv_orig = [min(mean(data_1))*1.1 max(mean(data_1))*0.7];%PSD_YLIM_ORG; %[-30,-10]; %zeros(1,2);
%             ylim_intv_fooof = [-1 max(mean(data_2))]*2; %[-2,6]; %zeros(1,2);
            
            figure('color','white','position',[100 300 600 240],'Renderer','Painters');
            j = 1;
            
            %{
            % -------------- original PSD -----------
            subplot(1,10,[j j+1])
            title('')
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
%             axsignif = highlight_CL(ax, fooof_freq, [], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('');
            ylabel('10*log(PSD)');
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
        %     title(['Cluster ',num2str(k)]);
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            subplot(1,10,[j+2 j+3])
            title('')
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
%             axsignif = highlight_CL(ax, fooof_freq, [], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'--','color','black');
            xlabel('');
            ylabel('');
            ax.YTick = [];
        %     set(ax,'LineWidth',1)
            set(ax,'FontName','Arial','FontSize',10)
        %     title(['Cluster ',num2str(k)]);
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');
            
            %}
            % -------------- Flattened PSD -----------
            subplot(1,5,[j j+1])
            title('');
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
            set(gca,'FontName','Arial','FontSize',10)
            ax = gca;       
%             axsignif = highlight_CL(ax, fooof_freq,[], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'-','color',[0.5 0.5 0.5]);
            ylabel('10*log(PSD)');
            xlabel('');
            set(gca,'FontName','Arial','FontSize',10)
        %     set(ax,'LineWidth',1)
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            subplot(1,5,[j+2 j+3])
            title('');
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
%             axsignif = highlight_CL(ax, fooof_freq,[], 'background', 'Frequency(Hz)');
            xlim([4 40]);
            plot([0 40],[0 0],'-','color',[0.5 0.5 0.5]);
            xlabel('');
            ylabel('');
            ax.YTick = [];
            set(gca,'FontName','Arial','FontSize',10)
        %     set(ax,'LineWidth',1)
            xline([4],'--'); xline([8],'--'); xline([13],'--'); xline([30],'--');

            
            subplot(1,5,[j+4]);hold on;
            color_stats = [223,101,176;66,146,198]/255;
            
            times = fooof_freq;
            for i = 1:2
                regions = pinter{k}{i}(:,2);
                in_a_region = 0;
                yl2(1) = 0-2*i;yl2(2) = -2-2*i;
                color2 = color_stats(i,:);
                for index=1:length(regions)
                    if regions(index) && ~in_a_region
                        tmpreg(1) = times(index);
                        in_a_region = 1;
                    end
                    if (~regions(index) || index == length(regions)) && in_a_region
                        tmpreg(2) = times(min(length(times), index));
                        in_a_region = 0;
                        tmph = patch([tmpreg(1) tmpreg(2) tmpreg(2) tmpreg(1)], ...
                            [yl2(1) yl2(1) yl2(2) yl2(2)], color2,'FaceAlpha',.5); hold on;
                        set(tmph, 'edgecolor', color2,'EdgeAlpha',.5);                       
                    end
                end
    
            end            
            xlim([4 40]);
            xlabel('Frequency (Hz)');
            ylabel('');
            set(gca,'FontName','Arial','FontSize',10)
            ax = gca;
            ax.YTick = [];
            ax.XTick = [10 20 30 40];
            ylim([-8 0]);
            box on;
        %     set(ax,'LineWidth',1)
            camroll(-270)
            j = j + 5;

%             if g == 1
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'.pdf']),'ContentType','vector','Resolution',300)
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'.jpg']),'ContentType','vector','Resolution',300)
%             elseif g == 2
%                 exportgraphics(gcf,fullfile(save_dir ,['Cluster_PSD_',num2str(k),'_speed.pdf']),'ContentType','vector','Resolution',300)
%             end
        end
    end
end

