function [resTable] = detectSpontRm(cellID,protNum,meta,resTableName)

%calculate estimated resistance during RC check (no tau, or Cm), and store
%baseline data

%created by Anne Takesian, Grace Chavez, and Christine Liu 8/2/24

%% load data and stim
[pathname] = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);

[traces,~] = loadData(pathname);
[allStimTable,~] = loadStim(pathname);
% 
% %Load required filter:
% load 'Lowpass equiripple filter.mat'

%% define stim parameters
dx = 0.1; %define acqusition freq (how many msec does each data point has)

% define variables related to RC check
if isempty(allStimTable.RC)
    RCcheck = 0;
    pulseStart = 100;
    pulseStop = pulseStart + 200;
else 
    RCcheck =1;
    pulseStart = allStimTable.RC.delay(1);
    pulseStop = pulseStart + allStimTable.RC.pulseWidth(1);
end 
  

%% Initialize variables
traceNum = size(traces,2); %number of data traces
% initialize arrays to store data 
RinRs =nan(traceNum,1);
peak =nan(traceNum,1);
baseline = nan(traceNum,1);


%%  calculate Rin
for i = 1:traceNum
   %Compute Rs and Rm for each trace
    if RCcheck
        [RinRs(i), peak(i), baseline(i)] = calcRs_VC_simple(traces(:,i), dx, allStimTable.RC(i,:));
    else 
        rin(i) = 0;
         rs(i) = 0;
         cm(i) = 0;
       error(i) = 0;
    end 

end
%% Generate result table

table_str = {'RinRs', 'peak', 'baseline'};
alldata_array = [RinRs,peak,baseline];
%if there were real data, convert the cell array into a mat to be stored in spikes table
resTable = array2table(alldata_array,'VariableNames',table_str);


%% save analysis
% save locally in current data folder
protPath = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);
save([protPath resTableName '.mat'],'resTable');


