function EEG = reject_eye_ic(EEG)
    fileout = [EEG.filepath filesep [EEG.filename(1:end-4),'_rej_eye.set']];
    % remove eye  ICs (classification >.9)
    rmEye = EEG.etc.ic_classification.ICLabel.classifications(:,3)>.9;
    EEG.ect.ICArej.IClabelEye = sum(rmEye);
%     procInfo(si).IC_eye = EEG.ect.ICArej.IClabelEye;

    % visualize
%     nComp = size(EEG.icawinv,2);
%     procInfo(si).IC_num = nComp;
%     EEG.setname = file_out;
%     rejComp = strjoin(string(find(rmEye)),', ');
%     pop_topoplot(EEG, 0, [1:nComp],...
%         [ID, ', rejected: ', rejComp{:}],...
%         [ceil(sqrt(nComp)) ceil(sqrt(nComp))] ,0,...
%         'electrodes','on','iclabel','on');
%     set(gcf, 'Units','normalized','Position',[0 0 1 1]);
%     print('-dpng', fullfile(PATHOUTsi, [ID, reportNames{1}])); % save
%     close;

    % reject components
    EEG = pop_subcomp( EEG, find(rmEye), 0); %add muscle here for ERP calc!

    % update info
    EEG.etc.ICA.ICremain = size(EEG.icawinv,2);
%     procInfo(si).IC_kept = EEG.etc.ICA.ICremain;

    % save time series data 
    pop_saveset(EEG,fileout);
    
%     fprintf('Done with sub-%03d!\n', si);
end