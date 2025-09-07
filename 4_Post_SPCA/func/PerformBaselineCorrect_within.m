function PerformBaselineCorrect_within(allersp,alltimes,allfreqs,config)
        alpha = 0.05;
        warpingvalues = config.warpingvalues;
        colormap_ersp = config.colormap_ersp;
        title_keyword = config.title_keyword;
        YlimMax = config.YlimMax;
        
        clear allerspdata_meanSubj
        allerspdata_meanSubj = allersp;
        
        % reorganize allerspdata
        clear erspdata baseline pboot
        for j = 1: length(allerspdata_meanSubj)
            erspdata = allerspdata_meanSubj{j}(:,:,:);
            baseidx = find(alltimes>=warpingvalues(1) & alltimes<=warpingvalues(5));                
            baseline_allcomp = mean(erspdata(:,baseidx,:),2); % mean power for each person
            baseline = mean(baseline_allcomp,3);%mean power across participant
            cluster_allcomp_ersp{j,1} = allerspdata_meanSubj{j}(:,:,:)-repmat(baseline_allcomp,1,length(alltimes));% subtract baseline for each person
%             cluster_allcomp_ersp_mean{j,1} = mean(allerspdata_meanSubj{j}(:,:,:)-repmat(baseline,1,length(alltimes)),3);
            cluster_allcomp_ersp_crop{j,1} = cluster_allcomp_ersp{j,1}(:,baseidx);
        end
%         std_plottf(alltimes2(baseidx),allfreqs2, cluster_allcomp_ersp_crop, 'datatype','ersp', 'plotmode','normal',...
%             'titles',alltitles,'caxis',[-climMat_max climMat_max])
%         saveas(gcf,fullfile(outputdir,['Component_'  '_ERSP_CLUSTER' , ' Design', file_keyword,num2str(STUDY.currentdesign)]));

        %% Paper Figure  - significance masked ERSP 
        freqvalues = [4 YlimMax];
        freqidx = find(allfreqs >= freqvalues(1) & allfreqs <= freqvalues(2));
        climMat = [min(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all') max(mean(cluster_allcomp_ersp{4}(freqidx,baseidx,:),3),[],'all')];
        climMat_max = max(abs(climMat));
        climMat_max = config.climMat_max;
        
        figure('color','white','position',[200 200 500 150],'renderer','Painters');
        for j = 1:length(allersp)%1:length(allersp)
            if ~isnan(alpha)
                fprintf('running the permutation\n');
                curr_ersp_temp = cluster_allcomp_ersp{j,1}(freqidx,baseidx,:);% this is already sub baseline
                curr_ersp_temp_mean = mean(curr_ersp_temp,3);
                surro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
                for n = 1:2000
                    bootLatency = randi(size(curr_ersp_temp,2),[size(curr_ersp_temp,2),1]); %random time sample
                    bootFreq = 1:size(curr_ersp_temp,1);
                    bootIc = 1:size(curr_ersp_temp,3); 
                    tmpSurro = mean(curr_ersp_temp(bootFreq,bootLatency,bootIc),3);
                    surro(:,:,n) = tmpSurro; %save 2000 iterations of surrogates 
                end

                bootSurro = zeros(size(curr_ersp_temp,1),size(curr_ersp_temp,2),2000);
                for n = 1:2000
                    bootIdx  = randi(2000,[size(curr_ersp_temp,3),1]);
                    tmpSurro = mean(surro(:,:,bootIdx),3);
                    bootSurro(:,:,n) = tmpSurro;
                end

                pvalMap = stat_surrogate_pvals(bootSurro,curr_ersp_temp_mean,'both');
                pvalMap(pvalMap>1)=1; 
                [p_masked, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalMap,0.05,'pdep',1);

                % debri removal
                [labelMap,uniqueLabelNum] = bwlabeln(p_masked);
                tmpDisp = sort(labelMap(:),'descend');
                [occurrence,idx] = hist(tmpDisp,unique(tmpDisp));
                sortOccurrence = sort(occurrence,'descend');
                disp(num2str(sortOccurrence(2:end)));
                threshold = 300;
                threshOccurrence = occurrence;
                threshIdx = find(threshOccurrence<threshold);
                kMask = ismember(labelMap,idx(threshIdx));
                finalMask = p_masked-kMask;

                clust_ersp = curr_ersp_temp_mean; 
                clust_maskedersp = clust_ersp; 
                clust_maskedersp(~finalMask) = 0;
    %             curr_maskedersp(~p_masked) = 0; 

            else
                clust_ersp = mean(cluster_allcomp_ersp_mean{j},3);
                clust_maskedersp = clust_ersp;
            end   
            
            subplot(1,length(allersp),j)
            colormap(colormap_ersp);
            faceAlpha_mask = ones(size(clust_maskedersp))*0.7; %alpha range [0-1] where 0 is fully transparent and 1 = fully opaque
            faceAlpha_mask(clust_maskedersp ~= 0 ) = 0; %0 is significant? 1 is not? 
        %     contourf(alltimes3, allfreqs3, erspDiff(4).raw,200,...
        %             'linecolor','none');hold on;
            contourf(alltimes(baseidx), allfreqs(freqidx), clust_ersp,200,...
                       'linecolor','none');hold on;
            imagesc(alltimes(baseidx),allfreqs(freqidx),clust_maskedersp,'AlphaData',faceAlpha_mask);
            %- add vertical line
            vline([warpingvalues(2) warpingvalues(3) warpingvalues(4)],{'k--' ,'k--', 'k--', 'k--'});
            set(gca,'clim',[-climMat_max, climMat_max],'xlim',[warpingvalues(1) warpingvalues(end)],...
                'ydir','norm','ylim',[4 YlimMax],'yscale','log')
            if YlimMax == 50
                set(gca,'YTick',[4,8,13,30,50]); 
                set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',10);
            elseif YlimMax == 100
                set(gca,'YTick',[4.01,8,13,30,50,99.4843]); 
                set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',10);
            end
            if j == 1
                ylabel(sprintf('Frequency (Hz)'),'fontsize',10,'fontweight','bold');
            else
                ylabel(sprintf(''),'fontsize',12,'fontweight','bold');
            end
            xlabel('','Fontsize',12);
%             title(strcat({'Cluster '},num2str(k),' ',terrain_keyword{j}));
            title(title_keyword{j});
%             set(gca,'xtick',warpingvalues,'xticklabel',{'RFS','LFO','LFS','RFO','RFS'});
            set(gca,'xtick',warpingvalues,'xticklabel',{'','','','',''});
            xtickangle(45)
            ax = gca; ax.XAxis.FontSize = 8;
        end
        hp4 = get(subplot(1,length(allersp),4),'Position');
        c = colorbar('Position',[hp4(1)+hp4(3)+0.01  hp4(2)-0.0085 0.008  hp4(4)-0.071]);
        c.Limits = [-climMat_max, climMat_max];
        hL = ylabel(c,{'\Delta Power(dB)'},'fontweight',...
            'bold','FontName','Arial','FontSize',8);
        set(hL,'Rotation',90);
        hL.Position(1) = 5;
    %             hL.Position(2) = hL.Position(2)+0.025;
        hL.Position(2) = hL.Position(2);
        c.Position(2) = .10;
        c.Position(4) = .5;    
%         colorbar
            
%         saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_within_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.fig']));
%         saveas(gcf,fullfile(output_folder,'allerspdata',['allerspdata_within_eeglab_',num2str(k),'_stats_',num2str(YlimMax),file_keyword,'.pdf']));
%         saveas(gcf,config.save_filepath);
end