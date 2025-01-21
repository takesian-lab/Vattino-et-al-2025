function [resTableName,resTable]= stim_updateSave(cellID,protNum,resTableName,resultTableDraft,traces,stimTraces,meta,setup,extraVars)
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
    if resultTableDraft.detect(i) == 0
        resultTableDraft(i,6:end) = emptyData([6:size(resultTableDraft,2)]);
    end 
end

%% filter data and traes based on Rs data in resultTableDraft --- 
% this is different step than replacing fields with NaN if there's a
% failure -- we are simply deleting the epochs that didn't meet the
% criteria for Rs requirements (which can be either pre-defined in the ops files 
% or defined here). [added 5/16/23]
if isfield(setup, 'RsThreshold')
    toKeepIndex = resultTableDraft.rs <= setup.RsThreshold;
else 
    toKeepIndex  = resultTableDraft.rs <= 100;
end 
% update resultTableDraft
resultTableDraft = resultTableDraft(toKeepIndex,:);

% update traces
traces = traces(:,toKeepIndex);

% update stimTraces
stimAOs = fieldnames(stimTraces);
for s = 1:length(stimAOs)
    if ~isempty(stimTraces.(stimAOs{s}))
       stimTraces.(stimAOs{s}) = stimTraces.(stimAOs{s})(toKeepIndex,:);
    end 
end 
%% Initialize variables
% this is different step than replacing fields with NaN if there's a
% failure -- we are simply deleting the epochs that didn't meet the
% criteria for Rs requirements (which can be either pre-defined in the ops files 
% or defined here). [added 5/16/23]
if isfield(setup, 'RsThreshold')
    toKeepIndex = resultTableDraft.rs <= setup.RsThreshold;
else 
    toKeepIndex  = resultTableDraft.rs <= 150;
end 
% update resultTableDraft
resultTableDraft = resultTableDraft(toKeepIndex,:);

traceNum = size(resultTableDraft,1);
dx = setup.dx; % sampling rate

preStimBaseline = resultTableDraft.preStimBaseline;
detect = resultTableDraft.detect;
onsetTime = resultTableDraft.onsetTime;

peakTimefromStim = resultTableDraft.peakTimefromStim; % or resultTableDraft.peakTimefromStim; depending on if output was already in ms
peakAmplitudeBaseline = abs(resultTableDraft.peakAmplitudeBaseline);
riseTime = resultTableDraft.riseTime;
totalCharge = abs(resultTableDraft.totalCharge);
fallTimefromPeak = resultTableDraft.fallTimefromPeak;
tauFall = resultTableDraft.tauFall;
rsquareMonoFit = resultTableDraft.rsquareMonoFit;
%% Tables generation:
successEvents = {sum(detect)};
totalTraceNum = {traceNum};

% write in extra vars
Var1 = extraVars.var1; var1Name = extraVars.var1Name;
Var2 = extraVars.var2; var2Name = extraVars.var2Name;
Var3 = extraVars.var3; var3Name = extraVars.var3Name;

Rin = {resultTableDraft.rin(~isnan(resultTableDraft.rin))};
Rs = {resultTableDraft.rs(~isnan(resultTableDraft.rs))};
Cm = {resultTableDraft.cm(~isnan(resultTableDraft.cm))};

Onset_time = {onsetTime(~isnan(onsetTime))};
Latency = {peakTimefromStim(~isnan(peakTimefromStim))};
Amplitude = {peakAmplitudeBaseline(~isnan(peakAmplitudeBaseline))};
RiseTime = {riseTime(~isnan(riseTime))};
Charge = {totalCharge(~isnan(totalCharge))};
FallTime = {fallTimefromPeak(~isnan(fallTimefromPeak))};
Tau_FallTime = {tauFall(~isnan(tauFall))};
Rsq_FallTime = {rsquareMonoFit(~isnan(rsquareMonoFit))};
resTable = table(protNum,successEvents,totalTraceNum,Var1,Var2,Var3,...
    Rin,Rs,Cm,Onset_time,Latency,Amplitude,RiseTime,Charge,FallTime,Tau_FallTime,Rsq_FallTime);
%replace all the empty cells with NaN
if resTable.successEvents{1} == 0
    resTable.Onset_time{1} = NaN;
    resTable.Latency{1} = NaN;
    resTable.Amplitude{1} = NaN;
    resTable.RiseTime{1} = NaN;
    resTable.Charge{1} = NaN;
    resTable.FallTime{1} = NaN;
    resTable.Tau_FallTime{1} = NaN;
    resTable.Rsq_FallTime{1} = NaN;
end
resTable.Traces = {traces};
resTable.stimTraces = {stimTraces};
resTable.Setup = {setup};

%% change the names of var1 and var2
if ~isempty(var1Name)
    var1Index =  find(contains(resTable.Properties.VariableNames,'Var1'));
    resTable.Properties.VariableNames{var1Index} = var1Name;
end 

if ~isempty(var2Name)
    var2Index =  find(contains(resTable.Properties.VariableNames,'Var2'));
    resTable.Properties.VariableNames{var2Index} = var2Name;
end

if ~isempty(var3Name)
    var3Index =  find(contains(resTable.Properties.VariableNames,'Var3'));
    resTable.Properties.VariableNames{var3Index} = var3Name;
end

%% save tables
% save this result table into current data folder
protPath = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);

save([protPath resTableName '.mat'],'resTable');


end