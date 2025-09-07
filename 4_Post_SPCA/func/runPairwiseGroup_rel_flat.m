function  [erspDiff_wind] = ...
    runPairwiseGroup_rel_flat(allersp,alltimes,allfreqs,config)   
    %% Compute the common baseline 
    clear baseline_allcomp baseline_allcond_mean baseline_allcond_median baseline_allcond_time baseline_allcond    
    fcn = @erspStats;
    fcn2 = @erspStats_V2;

    alpha               = 0.05;
    warpingvalues       = config.warpingvalues;
    colormap_ersp       = config.colormap_ersp;
    title_keyword       = config.title_keyword;     
    STUDY               = config.STUDY;
    runPairwiseCond     = config.runPairwiseCond;
    runPairwiseGroup    = config.runPairwiseGroup;
    k                   = config.k;   
    YlimMax             = config.YlimMax;
    climMat_max = [1];
    
    freqvalues = [4 YlimMax];
    baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5)); 
    freqidx = find(allfreqs>=freqvalues(1) & allfreqs<=freqvalues(2));
    
    
    % loop through each people, common baseline = mean(ERSP across
    % condition) for each person and then subtract for each person
    % allfreqs=200, alltimes=159
    for group = 1:size(allersp,2)
         clear allersp_nolog commen_ersp
         for n = 1:size(allersp{1,group},3)
            % convert to log
            for j = 1:4
                allersp_nolog(:,:,j,n) = 10.^(allersp{j,group}(:,:,n)/20);
            end
         end
%          common_ersp = mean(mean(mean(allersp_nolog(:,baseidx,:,:),3),4),2);
         common_ersp = mean(mean(allersp_nolog(:,baseidx,:,:),3),2);
         ref_ersp = mean(allersp_nolog(:,baseidx,1,n),2);
         for n = 1:size(allersp{1,group},3)
             for j = 1:4
                allersp_nolog_subBase{j,group}(:,:,n) = allersp_nolog(:,:,j,n)./repmat(common_ersp(:,:,n),1,size(allersp_nolog,2));
                allersp_log_subBase{j,group}(:,:,n) = 20*log10(allersp_nolog_subBase{j,group}(:,:,n)); 
             end
         end
    end
    
    allerspdata_meanSubj = allersp_log_subBase;
    
    fprintf('=== pairwise comparison === \n');
    f1 = figure('color','white','position',[200 200 700 200],'renderer','Painters');        
    for j = 2:length(allerspdata_meanSubj)
%             curr_ersp = allerspdata_meanSubj{j,1};
%             ref_ersp = allerspdata_meanSubj{j,2};
%             %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp ref_ersp});
%             [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp ref_ersp});
%             [erspDiff(j).raw] = [mean(curr_ersp,3)-mean(ref_ersp,3)];%young-old
%             [erspDiff(j).masked] = [erspDiff(j).raw.*pgroup{1,1}];
%             [erspDiff(j).pcond] = pgroup{1,1};
% 
        STUDY.etc.statistics.groupstats = 'on';
        STUDY.etc.statistics.condstats = 'off';
        %             STUDY.design.variable(1) = [];
        curr_ersp_wind = allerspdata_meanSubj{j,1}(freqidx,baseidx,:)-allerspdata_meanSubj{1,1}(freqidx,baseidx,:);%young
        ref_ersp_wind = allerspdata_meanSubj{j,2}(freqidx,baseidx,:)-allerspdata_meanSubj{1,2}(freqidx,baseidx,:);%old
        %[pcond, pgroup, pinter] = erspStats(STUDY,{curr_ersp;ref_ersp});
        [pcond, pgroup, pinter] = feval(fcn, STUDY,{curr_ersp_wind ref_ersp_wind});
        [erspDiff_wind(j).raw] = [mean(curr_ersp_wind,3)-mean(ref_ersp_wind,3)];
        [erspDiff_wind(j).masked] = [erspDiff_wind(j).raw.*pgroup{1,1}];
        [erspDiff_wind(j).pcond] = pgroup{1,1};

        subplot(1,length(allerspdata_meanSubj),j)
        colormap(colormap_ersp);
        faceAlpha_mask = ones(size(erspDiff_wind(j).pcond))*0.7; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
        faceAlpha_mask(erspDiff_wind(j).pcond == 1) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
        contourf(alltimes(baseidx), allfreqs(freqidx), erspDiff_wind(j).raw,200,...
                   'linecolor','none');hold on;
        imagesc(alltimes(baseidx), allfreqs(freqidx),erspDiff_wind(j).masked,'AlphaData',faceAlpha_mask);
        %- add vertical line
        vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
        set(gca,'xlim',[warpingvalues(1) warpingvalues(end)],...
            'ydir','norm','ylim',[allfreqs(1) YlimMax],'yscale','log')
        if YlimMax == 50
            set(gca,'YTick',[4,8,13,30,50]); 
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
        elseif YlimMax == 100
            set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
        end    
        if j == 1
            ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
        else
            ylabel(sprintf(''),'fontsize',10,'fontweight','bold');
        end
        xlabel('','Fontsize',10);
        title('');
        set(gca,'xtick',warpingvalues,'xticklabel',{'RFS','LFO','LFS','RFO','RFS'});
        xtickangle(45)
        ax = gca; ax.XAxis.FontSize = 8;
        ylim([4 YlimMax]);
        caxis([-climMat_max, climMat_max]);
    end
    hp4 = get(subplot(1,length(allerspdata_meanSubj),4),'Position');
    c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.012  hp4(4)]);
    c.Limits = [-climMat_max, climMat_max];
    hL = ylabel(c,[{'\Delta Power'};{'(dB)'}],'fontweight',...
        'bold','FontName','Arial','FontSize',9);
    set(hL,'Rotation',90);
    hL.Position(1) = hL.Position(1)-0.3;
    hL.Position(2) = hL.Position(2);

    exportgraphics(gcf,config.save_filepath);
end
