function plotERSP_parfor_crunch(tmpSTUDY,tmpSTUDY_commonbase,ALLEEG,myErspParams,clusternum,outputdir)
% Written by Chang - 2023-4-15 to run plotERSP with parallel processing
% output are saved as mat
% Update by Chang - 2024-01-24 for CRUNCH

%}
        % ----------------------------------------
        % In the std_erspplot_customParams, it also had the
        % std_readata_customParams. 
        disp('Gathering esrp after baseline correction using warp(1) to warp(5)')
        tic
        disp(tmpSTUDY.etc.statistics)
        disp(myErspParams);
        [tmpSTUDY,allerspdata2, alltimes2, allfreqs2, pgroup2,pcon2,pinter2] = std_erspplot_customParams(tmpSTUDY, ALLEEG, myErspParams, 'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_subBase.mat']), 'allerspdata2', 'alltimes2', 'allfreqs2', 'pgroup2', 'pcon2', 'pinter2');
        saveas(gcf,fullfile(outputdir,['Design_All_Comps_ERSP_subBase.fig']))
% 
%             toc

        % ---------------------------------------- Commonbase!!!
        disp('Gathering esrp after baseline correction using warp(1) to warp(5) and apply commonbase')
        tic
        disp(tmpSTUDY_commonbase.etc.statistics)
        disp(tmpSTUDY_commonbase.etc.erspparams)
        [tmpSTUDY_commonbase_out,allerspdata3, alltimes3, allfreqs3, pgroup3,pcon3,pinter3] = std_erspplot_customParams(tmpSTUDY_commonbase, ALLEEG, myErspParams, 'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_subBase_commonBase.mat']), 'allerspdata3', 'alltimes3', 'allfreqs3', 'pgroup3', 'pcon3', 'pinter3');
        saveas(gcf,fullfile(outputdir,['Design_All_Comps_ERSP_subBase_commonBase.fig']))
% 
                        % ---------------------------------------- Commonbase!!!
%}
%{
        disp('Gathering esrp after fulltrial correction and baseline correction using warp(1) to warp(5) and apply commonbase')
        tic
        disp(tmpSTUDY_commonbase.etc.statistics)
        disp(tmpSTUDY_commonbase.etc.erspparams)
        % if trialbase = 'on', no commonbaseline subtraction is
        % performed in newtimefbaseln
        [tmpSTUDY_commonbase_out,allerspdata4, alltimes4, allfreqs4, pgroup4,pcon4,pinter4] = std_erspplot_customParams(tmpSTUDY_commonbase, ALLEEG, myErspParams_trialbasefull,...
            'clusters', clusternum,...
            'freqrange', [0 200]); 
        save(fullfile(outputdir,['readESRP_', num2str(STUDY.currentdesign),file_keyword,'_subBase_fullTrialbase_commonBase.mat']), 'allerspdata4', 'alltimes4', 'allfreqs4', 'pgroup4', 'pcon4', 'pinter4');
        saveas(gcf,fullfile(outputdir,['Design',num2str(STUDY.currentdesign),file_keyword, '_All_Comps_ERSP_subBase_fullTrialbase_commonBase.fig']))
        %}
        % 
%}
end