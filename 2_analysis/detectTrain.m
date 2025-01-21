function [resTableName,resTable] = detectTrain(cellID,protNum,meta,polarity,Setup,stimType,drug)

%Dongqin May 2019; edited by Christine J. Liu Jan 2021
%This program is designed to analyze synaptic currents evoked by train
%stimulation with both 10 Hz and 20Hz

%% load data and stim
[pathname] = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);

[traces,~] = loadData(pathname);
[allStimTable,stimTraces] = loadStim(pathname);

%Load required filter:
load 'Lowpass equiripple filter.mat'

%% define stim parameters

% define variables from the Setup struct
startSearch = Setup.startSearch; %when to stop searching for peak after stim
stopSearch = Setup.stopSearch; %time before artifact to average baseline
dx = Setup.dx;

AOfields = fieldnames(allStimTable);
aoInd = contains(AOfields,stimType,'IgnoreCase',true);
if sum(aoInd) == 1
    StimParams_all = allStimTable.(AOfields{aoInd});
    StimTraces_all = stimTraces.(AOfields{aoInd});
elseif sum(aoInd) == 0
    disp('no AO fields match your stimType entry. Please recheck.');
    disp(AOfields);
else 
    disp('the following AO fields matches your stimType entry:')
    disp(AOfields(aoInd))
    newStimType = input('Which one would you like to analyze now?...');
    pause;
    aoIndNew = contains(AOfields,newStimType,'IgnoreCase',true);
    StimParams_all = allStimTable.(AOfields{aoIndNew});
    StimTraces_all = stimTraces.(AOfields{aoIndNew});
end 

%make sure that all the stims are the same for each acq within this protocal
if size(unique(StimParams_all),1)==1 %if there's only 1 unique row for this stim table
    StimParams = unique(StimParams_all);
    StimTracesTemp = find(StimTraces_all(1,:)>0);
else 
    disp('There are more than one type in this protocal - not compatible with this analysis type.');
end 

%Parameters for the stimulation protocol (in ms):
delay = StimParams.delay;
numStim = StimParams.numPulses;
isi = StimParams.isi;
pulserate = 1/isi * 1000; %in Hz

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
traceNum = size(traces,2);

% Define experiment-specific parameters
sampling = 1:size(traces,1); %sampling points: 10000 points per second
msTime = 1/dx; %for converting points to ms --> 1/dx is 10 points = 1 ms
time = sampling/msTime; %in ms, given the ms converter

baselineStart = 1;
baselineEnd = pulseStart-1;
pulseLength = pulseStop-pulseStart;

pulsePoint = 1/dx*(delay:isi:delay+((numStim-1)*isi)); %defines the artifact times taking into account the freq
vectorofMax = zeros(1,numStim);
indexofMax = zeros(1,numStim);
traceNum = size(traces,2);
globalPeak = NaN(traceNum,numStim);
localPeak = NaN(traceNum,numStim);
rin = NaN(1,traceNum);
rs = NaN(1,traceNum);
cm = NaN(1,traceNum);
preStimBaseline = mean(traces(pulsePoint(1)-100*msTime:pulsePoint(1),1));

%% plot trains and calculate
figure;
for i = 1:traceNum
   %Compute Rs and Rm for each trace
    if RCcheck
        [rin(i), rs(i), cm(i), error(i)] = calcRs_VC(traces(:,i), dx, allStimTable.RC(i,:));
    else 
        rin(i) = 0;
         rs(i) = 0;
         cm(i) = 0;
       error(i) = 0;
    end 
    
    globalBaseline(i) = mean(traces(pulsePoint(1)-100*msTime:pulsePoint(1),i)); %uses last 10 ms before 1st pulse to compute baseline
    for j = 1:numStim % number of pulses
        if polarity == 1
            [vectorofMax(j), indexofMax(j)] = max(traces(pulsePoint(j)+startSearch:pulsePoint(j)+stopSearch,i)); % +50 to avoid artifacts;
        elseif polarity == -1
            [vectorofMax(j), indexofMax(j)] = min(traces(pulsePoint(j)+startSearch:pulsePoint(j)+stopSearch,i));
            
        end
        indexofMax(j) = indexofMax(j)+ pulsePoint(j) + startSearch-1;
        
        localBaseline(j) = mean(traces(pulsePoint(j)-msTime:pulsePoint(j),i)); %1/dx is 10 points = 1 ms
        
        localPeak(i,j) = abs(vectorofMax(j)-localBaseline(j));
        
    end
    
    % indexofMax is the time where the max was found
    % vectorofMax is the amplitude of that max
    
    globalPeak(i,:) = abs(vectorofMax-globalBaseline(i)); %calculates the amplitude as the difference between the max and the baseline before 1st pulse
    %localPeak(i,:) = abs(vectorofMax(j)-localBaseline(j));
    
    subplot1 = subplot(traceNum,1,i);
    plot(time,traces(:,i));
    hold on
    plot(indexofMax/msTime,vectorofMax,'ro');
    
    if polarity == 1
        ymax = max(vectorofMax)+40;
        ylim(subplot1, [preStimBaseline-50 ymax]);
    else
        ymin = min(vectorofMax)-40;
        ylim(subplot1, [ymin preStimBaseline+200]);
    end
    xlim(subplot1, [(pulsePoint(1)-1500)/msTime, (pulsePoint(end)+10000)/msTime]);
    axis off
       
end

% figname=[prot_result,'_figure','.fig'];
% saveas(gcf, fullfile(pathname,figname));

suppressionindex = mean(mean(globalPeak(end-1:end,:)))/mean(globalPeak(1,:)); %do the same for localPeak?
suppressionindex_all = suppressionindex*ones(1,traceNum);

if isfield(Setup, 'RsThreshold')
    toKeepIndex = rs <= Setup.RsThreshold;
else 
    toKeepIndex  = rs <= 100;
end 

% update traces
traces = traces(:,toKeepIndex);

% update stimTraces
stimAOs = fieldnames(stimTraces);
for s = 1:length(stimAOs)
    if ~isempty(stimTraces.(stimAOs{s}))
       stimTraces.(stimAOs{s}) = stimTraces.(stimAOs{s})(toKeepIndex,:);
    end 
end 
%% Generate result table
suppressionIndex = {suppressionindex_all(toKeepIndex)};
Rin = {rin(toKeepIndex)};
Rs = {rs(toKeepIndex)};
Cm = {cm(toKeepIndex)};
Drug = {drug};
globalPeak_toKeep = globalPeak(toKeepIndex,:);
localPeak_toKeep = localPeak(toKeepIndex,:);
GlobalPeakAmp = {array2table(globalPeak_toKeep)};%'VariableNames',{'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10'})};
LocalPeakAmp = {array2table(localPeak_toKeep)};


resTable = table(protNum,polarity, Drug, pulserate, suppressionIndex,Rin,Rs,Cm,GlobalPeakAmp,LocalPeakAmp);
resTable.Traces = {traces};
resTable.stimTraces = {stimTraces};
resTable.Setup = {Setup};

%% save analysis
% save locally in current data folder
resTableName = 'trainResults';
protPath = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);
save([protPath resTableName '.mat'],'resTable');


