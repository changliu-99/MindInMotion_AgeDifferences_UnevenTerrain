function [allersp, alltimes, allfreqs] = compute_ersp(STUDY,EEG)

    icatimf_f = [EEG.filepath filesep sprintf('%s.icatimef',EEG.subject)];

%     TIMEWARP_NTIMES = floor(EEG(1).srate/pi); % conservative nyquist frequency. making this too big can cause overlap between gait cyles

    filename = dir(fullfile(EEG.filepath,'*.icatimef'));
    load([filename.folder filesep filename.name],'-mat','parameters');
    
    if ~isempty(find(strcmp(parameters,'timewarpms')))
        warpingvalues = round(parameters{find(strcmp(parameters,'timewarpms'))+1});
        baseline_range = [warpingvalues(1) warpingvalues(end)];    
        myErspParams = parameters;
    %   myErspParams{find(strcmp(parameters,'baseline'))+1} = baseline_range; %'baseline';
        myErspParams{length(parameters)+1} = 'trialbase';
        myErspParams{length(parameters)+2} = 'off';%I don't think it works for trial 'on' condition somehow
    else
        baseline_range = [];
        myErspParams = parameters;
    %   myErspParams{find(strcmp(parameters,'baseline'))+1} = baseline_range; %'baseline';
        myErspParams{length(parameters)+1} = 'trialbase';
        myErspParams{length(parameters)+2} = 'off';%I don't think it works for trial 'on' condition somehow
        myErspParams{length(parameters)+3} = 'timewarpms';
        myErspParams{length(parameters)+4} =  'NaN';
    end
    cellArray= {myErspParams{1,2:2:length(myErspParams)}};
    fields = {myErspParams{1,1:2:length(myErspParams)-1}};
    myErspParams = cell2struct(cellArray, fields,2);


    [STUDY, allersp, alltimes, allfreqs, events, paramsersp] = std_readdata_customParams(STUDY, EEG, myErspParams, ...
        'datatype', 'ersp', ...
        'components', 1:size(EEG.icawinv,2),'clusters',1, 'singletrials', 'off', 'subbaseline', 'off', ...
        'timerange', baseline_range, 'freqrange', [4,60], ...
        'design', 1); %, 'concatenate', params.concatenate
     
    save([EEG.filepath filesep 'allersp.mat'],'allersp','alltimes');
end