function [resTableName,resTable]= minStim_updateSave(cellID,protNum,stimThreshold,resultTableDraft,traces,stimTraces,meta,monoSetup)
% updates and saves the analyzed optoStim data into projData defined in
% meta struct and saves the result table in the corresponding data folder.
% Should be used after running "detectStim" function on the same
% protocal.

% Updated on 2023/5/16 by CJL
% TODO: try to reduce the number of inputs to this function

%% Make sure that when monosynaptic == 0, all values are NaN!
emptyData = cell(1,size(resultTableDraft,2)); %to create empty data to call for replacing valuces in the cell;
emptyData(:)={NaN};
for i = 1:size(resultTableDraft,1)
    if resultTableDraft.detect_mono(i) == 0
        resultTableDraft(i,6:end) = emptyData([6:size(resultTableDraft,2)]);
    elseif resultTableDraft.detect_disy(i) == 0
        resultTableDraft(i,10:end) = emptyData([10:size(resultTableDraft,2)]);
    end 
end

%% filter data and traes based on Rs data in resultTableDraft --- 
% this is different step than replacing fields with NaN if there's a
% failure -- we are simply deleting the epochs that didn't meet the
% criteria for Rs requirements (which can be either pre-defined in the ops files 
% or defined here). [added 5/16/23] 
% - removed 7/19/23 - can be filtered out later during analysis 


if isfield(monoSetup, 'RsThreshold')
    toKeepIndex = resultTableDraft.rs <= monoSetup.RsThreshold;
else 
    toKeepIndex  = resultTableDraft.rs <= 100;
end 
% update resultTableDraft
resultTableDraft_updated = resultTableDraft;
placeholder = NaN(sum(~toKeepIndex),size(resultTableDraft,2));
resultTableDraft_updated(~toKeepIndex,:) = array2table(placeholder);
% 
% % update traces
% traces = traces(:,toKeepIndex);
% 
% % update stimTraces
% stimAOs = fieldnames(stimTraces);
% for s = 1:length(stimAOs)
%     if ~isempty(stimTraces.(stimAOs{s}))
%        stimTraces.(stimAOs{s}) = stimTraces.(stimAOs{s})(toKeepIndex,:);
%     end 
% end 
%% Initialize variables

traceNum = size(resultTableDraft,1);
dx = monoSetup.dx; % sampling rate

preStimBaseline = resultTableDraft.preStimBaseline;
detect_mono = resultTableDraft.detect_mono;


onsetTime = resultTableDraft_updated.onsetTime;
finalFirstPeakTime = resultTableDraft_updated.finalFirstPeakTime; % or resultTableDraft.peakTimefromStim; depending on if output was already in ms
finalFirstPeakAmp = abs(resultTableDraft_updated.finalFirstPeakAmp);
detect_disy = resultTableDraft_updated.detect_disy;
finalSecPeakTime = resultTableDraft_updated.finalSecPeakTime;
finalSecPeakAmp = resultTableDraft_updated.finalSecPeakAmp;

%% Tables generation:

minStimThreshold = {stimThreshold};
minTraces = {traces};
minstimTraces = {stimTraces};
minSetup = {monoSetup};
successMonos = {sum(detect_mono)};
totalTraceNum = {traceNum};

Rin = {resultTableDraft.rin};
Rs = {resultTableDraft.rs};
Cm = {resultTableDraft.cm};

% Onset_time = {onsetTime(~isnan(onsetTime))};
% monoLat = {finalFirstPeakTime(~isnan(finalFirstPeakTime))};
% monoAmp = {finalFirstPeakAmp(~isnan(finalFirstPeakAmp))};
% detect_disy = {detect_disy(~isnan(detect_disy))};
% diLat = {finalSecPeakTime(~isnan(finalSecPeakTime))};
% diAmp = {finalSecPeakAmp(~isnan(finalSecPeakAmp))};
detect_mono = {detect_mono};
Onset_time = {onsetTime};
monoLat = {finalFirstPeakTime};
monoAmp = {finalFirstPeakAmp};
detect_disy = {detect_disy};
diLat = {finalSecPeakTime};
diAmp = {finalSecPeakAmp};

resTable = table(protNum,minStimThreshold,minTraces,minstimTraces,minSetup,...
    successMonos,totalTraceNum,Rin,Rs,Cm,...
    detect_mono,Onset_time,monoLat,monoAmp,detect_disy,diLat,diAmp);
%replace all the empty cells with NaN
if resTable.successMonos{1} == 0
    resTable.Onset_time{1} = NaN;
    resTable.monoLat{1} = NaN;
    resTable.monoAmp{1} = NaN;
    resTable.detect_disy{1} = NaN;
    resTable.diLat{1} = NaN;
    resTable.diAmp{1} = NaN;
end

%% save tables
% save this result table into current data folder
resTableName = 'minResults';
protPath = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);
save([protPath resTableName '.mat'],'resTable');

end