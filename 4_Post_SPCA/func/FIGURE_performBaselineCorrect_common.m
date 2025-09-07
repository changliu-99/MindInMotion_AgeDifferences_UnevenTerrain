function FIGURE_performBaselineCorrect_common_tworow(allersp,alltimes,allfreqs,pcond,pgroup,pinter,config)

    alpha               = 0.05;
    warpingvalues       = config.warpingvalues;
    colormap_ersp       = config.colormap_ersp;
    STUDY               = config.STUDY;
    k                   = config.k;   
    YlimMax             = config.YlimMax;
    climMat_max         = config.climMat_max;
    title_keyword       = config.title_keyword;

    freqvalues = [4 YlimMax];
    baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5)); 
    freqidx = find(allfreqs>=freqvalues(1) & allfreqs<=freqvalues(2));

%%  Figure: tfplot with preset range
    figure('color','white','position',[200 200 750 600],'renderer','Painters');
    p = 1;
    colormap(colormap_ersp);
    for group = 1:2
        for j = 1:length(allersp)
            clust_ersp = mean(allersp{j,group},3);
            clust_maskedersp = clust_ersp;

            subplot(3,length(allersp)+1,p)
            tftopo(clust_maskedersp,alltimes(baseidx),allfreqs(freqidx),'limits',... 
                [warpingvalues(1) warpingvalues(end) nan nan -climMat_max climMat_max],...
                'vert',[warpingvalues(2:4)],'logfreq','native');
    %             ylim(log([3 50]));    
            ylim(log([4 YlimMax]))
            if YlimMax == 50
                set(gca,'YTick',log([4.01,8,13,30,50])); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
            elseif YlimMax == 100
                set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
            end    
            set(gca,'clim',[-climMat_max, climMat_max]);
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
                set(gca,'YTickLabel',{'','','','',''},'Fontsize',10);
            end
            xlabel('','Fontsize',10);
            set(gca,'XTick',warpingvalues,'XTicklabel',{'RFS','LFO','LFS','RFO','RFS'});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
            if group == 1
                title(title_keyword{j}); 
            else
                title('');
            end
            xlim([warpingvalues(1)-0.01 warpingvalues(end)+0.01]);
            colormap(colormap_ersp);

            % stat
            p = p + 1;
        end
        p = p + 1;
    end
    % --- add stats plot ------------
    subplot(3,length(allersp)+1,5) % add one subplot for stats
    tftopo(double(pcond{1}),alltimes(baseidx),allfreqs(freqidx),'limits',... 
        [warpingvalues(1) warpingvalues(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warpingvalues(2:4),'logfreq','native');
%             ylim(log([3 50]));
    colormap(colormap_ersp);
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'','','','',''},'Fontsize',10);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
    end    

    set(gca,'XTick',warpingvalues,'XTicklabel',{'','','','',''});
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 8;
    title('F Stats');  
%     title(strcat({' '},num2str(k)));
    ylabel('');xlabel('');xlim([warpingvalues(1)-0.01 warpingvalues(end)+0.01]);

%             cbar('vert',1:64,[-climMat_max, climMat_max])
%     hp4 = get(subplot(3,length(allersp)+1,5),'Position');
%     c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)-0.071]);
%     c.Limits = [-climMat_max, climMat_max];
%     hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
%         'bold','FontName','Arial','FontSize',9);
%     set(hL,'Rotation',0);
%     hL.Position(1) = hL.Position(1)+1.7;
%     hL.Position(2) = hL.Position(2)+0.025;
%     hL.Position(2) = .13;
%     c.Position(2) = .13;
%     c.Position(4) = .73;

     % --- add stats plot ------------
    subplot(3,length(allersp)+1,10) % add one subplot for stats
    tftopo(double(pcond{2}),alltimes(baseidx),allfreqs(freqidx),'limits',... 
        [warpingvalues(1) warpingvalues(end) nan nan  [-climMat_max, climMat_max]],...
        'vert',warpingvalues(2:4),'logfreq','native');
%             ylim(log([3 50]));
    ylim(log([4 YlimMax]))
    if YlimMax == 50
        set(gca,'YTick',log([4.01,8,13,30,50])); 
        set(gca,'YTickLabel',{'','','','',''},'Fontsize',10);
    elseif YlimMax == 100
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843])); 
        set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
    end    

    set(gca,'XTick',warpingvalues,'XTicklabel',{'RFS','LFO','LFS','RFO','RFS'});
    xtickangle(45)
    ax = gca; ax.XAxis.FontSize = 8;
    title('');  
%     title(strcat({'F stats Cluster '},num2str(k)));
    colormap(colormap_ersp);
    ylabel('');xlabel('');xlim([warpingvalues(1)-0.01 warpingvalues(end)+0.01]);

%             cbar('vert',1:64,[-climMat_max, climMat_max])
    hp4 = get(subplot(3,length(allersp)+1,10),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',90);
    hL.Position(1) = hL.Position(1)-0.3;
    hL.Position(2) = hL.Position(2);
%     hL.Position(2) = .13;
%     c.Position(2) = .13;
%     c.Position(4) = .73;
    
    for cond = 1:4
        curr_ersp_wind = allersp{cond,1}(:,:,:);%young
        ref_ersp_wind = allersp{cond,2}(:,:,:);%old
        [erspDiff_wind(cond).raw] = [mean(curr_ersp_wind,3)-mean(ref_ersp_wind,3)];
        [erspDiff_wind(cond).masked] = [erspDiff_wind(cond).raw.*pgroup{1,cond}];
        [erspDiff_wind(cond).pgroup] = pgroup{1,cond};
        
        subplot(3,length(allersp)+1,cond+10) % add one subplot for stats
        faceAlpha_mask = ones(size(erspDiff_wind(cond).pgroup))*0.6; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        faceAlpha_mask(erspDiff_wind(cond).pgroup == 1) = 0; %0 is significant? 1 is not? 
    %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
    %             'linecolor','none');hold on;
        contourf(alltimes(baseidx), allfreqs(freqidx), erspDiff_wind(cond).raw,200,...
                   'linecolor','none');hold on;
        imagesc(alltimes(baseidx),allfreqs(freqidx),erspDiff_wind(cond).masked,'AlphaData',faceAlpha_mask,[-climMat_max, climMat_max]);
        vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
        set(gca,'xlim',[warpingvalues(1) warpingvalues(end)],...
            'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
        caxis(c.Limits)
        if YlimMax == 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
        elseif YlimMax == 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
        end
        xlabel('','Fontsize',8);
        set(gca,'XTick',warpingvalues,'XTicklabel',{'RFS','LFO','LFS','RFO','RFS'});
        xtickangle(45)
        if cond == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
            set(gca,'YTickLabel',{'','','','',''},'Fontsize',10);
        end
        ax = gca; ax.XAxis.FontSize = 8;ylim([4 50]);
        colormap(colormap_ersp);
    end
    %     if ~isnan(alpha)
    %         saveas(gcf,fullfile(output_folder,'allerspdata3',['allerspdata3_eeglab_2_',num2str(k),'_stats_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    %         saveas(gcf,fullfile(output_folder,'allerspdata3',['allerspdata3_eeglab_2_',num2str(k),'_stats_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    %     else
    %         saveas(gcf,fullfile(output_folder,'allerspdata3',['allerspdata3_eeglab_2_',num2str(k),'_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    %         saveas(gcf,fullfile(output_folder,'allerspdata3',['allerspdata3_eeglab_2_',num2str(k),'_notfull_',num2str(YlimMax),file_keyword,'.fig']));
    %     end
    exportgraphics(gcf,config.save_filepath);
    
end