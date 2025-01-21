function [resTableName,resTable] = detectStim_CC(cellID,protNum,meta,polarity,Setup,stimType,extraVars,resTableName)

% Function to analyze single peak (ex. minimal stim protocal for thalamic stim project, or single opto stim analysis - as long as
% we don't look at kinetics)
% This script will go through each trace, and detect if there are any
% single or double peaks, and will output a resultsTableDraft, which can be
% changed and saved later.
% Inputs:
% 1. cellID - use "cell1" format (no space)
% 2. protNum - use number format
% 3. meta - setting structure that has all the info about the recording day/cell
% 4. polarity - value indicating if IPSCs or EPSCs
% 5. optoSetup - structure includes predefined changable variables
% Outputs:
% 1. resultTableDraft: output a table for to be later edited per protocal.
% 2. traces: the filtered traces to be later accessed and saved

%Created by Lucas & Christine 2020;
%Lucas turned into function in Nov 2022
%Edited by Christine + Lucas for github in Nov 2022
%Updated to simplify loading data and stim & updated calcRs_CC integration by Christine on 3/10/23

%% load data and stim
[pathname] = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);

[traces,~] = loadData(pathname);
[allStimTable,stimTraces] = loadStim(pathname);

%Load required filter:
load 'Lowpass equiripple filter.mat'

%% define stim parameters
% make sure that all the stims are the same for each acq within this protocal
AOfields = fieldnames(allStimTable);
aoInd = contains(AOfields,stimType,'IgnoreCase',true);
if sum(aoInd) == 1
    StimParams = allStimTable.(AOfields{aoInd});
elseif sum(aoInd) == 0
    disp('no AO fields match your stimType entry. Please recheck.');
    disp(AOfields);
else
    disp('the following AO fields matches your stimType entry:')
    disp(AOfields(aoInd))
    newStimType = input('Which one would you like to analyze now?...');
    aoIndNew = contains(AOfields,newStimType,'IgnoreCase',true);
    StimParams = allStimTable.(AOfields{aoIndNew});
end

% check if all the StimParams are the same for each acq, if they are, just
% store them as one single row in the stim table
if size(unique(StimParams),1) == 1
    StimParams = unique(StimParams);
end

%Parameters for the stimulation protocol (in ms):
stimStart = StimParams.delay;
searchStart = stimStart + Setup.searchAfterStim; %when to stop searching for peak after stim
searchStop = searchStart + Setup.searchWindow; %when to stop searching for peak after stim
preStim = stimStart - Setup.timePreStim; %time before artifact to average baseline

if isempty(allStimTable.RC)
    RCcheck = 0;
else
    RCcheck =1;
    pulseStart = allStimTable.RC.delay(1);
    pulseStop = pulseStart + allStimTable.RC.pulseWidth(1);
end


%% ---- load from Setup for analysis specific variables
%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
dx = Setup.dx;
% -- for finding multiple peask 
minPeakProm = Setup.minPeakProm;
minPeakDist = Setup.minPeakDist;

numPeaks = Setup.numPeaks;

peakAmplitudeRef = Setup.peakAmplitudeRef;


%% Initialize variables:
traceNum = size(traces,2);

baseline = NaN(1,traceNum);
preStimBaseline = NaN(1,traceNum);

%Turn time into points using sampling rate dx:
datastimStart = stimStart/dx;
datasearchStart = searchStart/dx;
datasearchStop = searchStop/dx;
datapreStim = preStim/dx;


onsetTime_ref_peak1 = Setup.onsetTime_ref_peak1 + datastimStart;
onsetTime_ref_peak2 = Setup.onsetTime_ref_peak2 + datastimStart;

%----

sampling = 1:size(traces,1); %sampling points
time = sampling*dx;

%% Analyze traces:

figure('units','normalized'); %creates fig that will be filled with subplot
% initialize variables
detectAP = zeros(1,traceNum); %logical array that identifies traces with AP

% Double arrays to store the results for mono (Peak1) and di (Peak2)
% synaptic events (if the peaks meet the critria for mono and disynaptic
% events
Peak1_OnsetLat = nan(1,traceNum);
Peak1_OnsetAmp = nan(1,traceNum);
Peak1_Amp = nan(1,traceNum);
Peak1_Lat = nan(1,traceNum);

Peak2_OnsetLat = nan(1,traceNum); %msec after stim onset
Peak2_OnsetAmp = nan(1,traceNum); %real amp NOT subtracting from baseline 
Peak2_Amp = nan(1,traceNum); %msec after stim onset
Peak2_Lat = nan(1,traceNum); %amp subtracting from prestim baseline 

% store all the detected peaks in cell arrays for each cell
all_PeakLat = cell(1,traceNum); %cell array to store all the cells 
all_PeakAmp = cell(1,traceNum);

% store max slope for the first peak (rise) 
maxLeftSlope = nan(1,traceNum);
% store max peak amp from baseline
maxPeakAmp = nan(1,traceNum);

% store AUC for the trace
PSP_auc = nan(1,traceNum);

%% Trace by trace analysis:
for i=1:traceNum %specifies that we refer to 2nd col: 50 (# of traces) and not to the 1st (number of time points)
    %smooth and filter each trace;
    data = smooth(filter(Hd1,traces(:,i))); %Hd1 is custom made filter, loaded at the beggining
    
    %Compute Rs and Rm for each trace
    if RCcheck
        [rin(i), rs(i), cm(i), error(i)] = calcRs_CC(traces(:,i), dx, allStimTable.RC(i,:));
    else
        [rin(i), rs(i), cm(i), error(i)]  = deal(NaN);
    end
    
    
    %% find peaks in data (peaks in EPSP or toughs in IPSP)
    
    %compute baselines for each trace
%     baseline(i) = mean(data(1:pulseStart-1)); % baseline before Rm step
    preStimBaseline(i)= mean(data(datapreStim:datastimStart-1)); %baseline 50 ms before stimulus
    
    %compute peak of synaptic current, taking into account if EPSP or IPSP:
    if polarity == -1
        data = -data;
    end
    
    % find multiple peaks: amplitude and location (in data points)
    [peakAmp, peakPoint] = findpeaks(data(datasearchStart:datasearchStop-1),...
        'MinPeakDistance',minPeakDist,'MinPeakProminence',minPeakProm);
    dataPeakPoint = peakPoint + datasearchStart;
    realPeakAmp = peakAmp - preStimBaseline(i);
   
        
    % first check for amplutides: get rid of any that didn't meet
    % peakAmplitudeRef in set up
    peakAmp(realPeakAmp<peakAmplitudeRef) = [];
    peakPoint(realPeakAmp<peakAmplitudeRef) = [];
    dataPeakPoint(realPeakAmp<peakAmplitudeRef) = [];
    realPeakAmp(realPeakAmp<peakAmplitudeRef) = [];
    
    % write all the onset amp and lat from this trace into the cell array
    % of all traces
    all_PeakLat{i} = peakPoint*dx; %in msec
    all_PeakAmp{i} = realPeakAmp; % amp in mV subtracting baseline
    if isempty(realPeakAmp)
        continue
    else 
    maxPeakAmp(i) = max(realPeakAmp);
    
    % if there are more than one peaks
    if length(peakPoint)>1
        % find the troughs between each of the peaks
        datatroughPoint = nan(length(peakPoint)-1,1); %establissh toughTime (should be one less than number of peaks)
        troughAmp = nan(length(peakPoint)-1,1); %establissh toughTime (should be one less than number of peaks)
        for t = 1:length(datatroughPoint)
            dataBW = data(dataPeakPoint(t):dataPeakPoint(t+1));
            troughAmp(t)=min(dataBW); %
            datatroughPoint(t) = dataPeakPoint(t) + find(dataBW == troughAmp(t)); % in data points since the start of trace
        end
    end % end of if length(peakTime) >1
    
    
    %% calculate onset for each peak in this trace
    % initialize onsetAmp and onsetTime vairables
    onsetAmp = nan(length(peakPoint),1);
    onsetPoint = nan(length(peakPoint),1); % number of datapoints after stim
    
    for n = 1:length(peakPoint)
        if n == 1 % first peak
            if peakAmp(n) > 0 % if this peak is an AP
                detectAP(i) = 1;
                [onsetAmp(n),onsetPoint(n)] = findOnsetAP(data,datasearchStart,dataPeakPoint(n),Setup.limitAPthresh);
            else
                [onsetAmp(n),onsetPoint(n)] = findOnsetTime(data,dataPeakPoint(n),stimStart./dx,Setup);
            end
        else 
            thisStimStart = datatroughPoint(n-1)- Setup.searchAfterStim./dx; %look for onset time after the trough right before this peak 
            [onsetAmp(n),onsetPoint(n)] = findOnsetTime(data,dataPeakPoint(n),thisStimStart,Setup);       
        end
    end %end looping through each
    
    
    %% find first peak max slope 
    % find the first derivative of the data between the start of data seach
    % and the first onset point (left slope)
    ix = datasearchStart:dataPeakPoint(1);
    maxLeftSlope(i) = max(diff(data(ix),1)); 
    
    %% AUC for all EPSPs
    
    % ---- comput fall time after the last peak     
    % from the last peak, find 3 consecutive points that are smaller than
    % 10% of last peak amp from baseline 
    lastdata = data(dataPeakPoint(end):end);
    tenPercentAmplitude = (preStimBaseline(i)+(0.1*realPeakAmp(end)));
    returnedDataInd = find(lastdata<tenPercentAmplitude); % indices for all the data points that are smaller than 10% amp
    if ~isempty(returnedDataInd)
        fallPoint = dataPeakPoint(end) + returnedDataInd(3); % find the third point and add the point to the last peak data point 
    else
        fallPoint = length(data); 
        disp('no synaptic current fall detected');
    end 
    
    %--- define the x and y for the trapz function 
    y = data-preStimBaseline(i); %center data along prestim baseline 
    t = time./1000;  % re-define time in seconds (for later intergral)
    
    % ---define period in data points that we are intergrating 
    p_ix = onsetPoint(1):fallPoint; 
    % TODO: make the onset detection more realizble so that we can use
    % first peak onset as the start to calculate integral 
    
    % ---find time-integral of voltage (mV*sec)
    PSP_auc(i) = trapz(t(p_ix),y(p_ix));
    
    %% check if the peak onset, time, and amp qualifies for the desired number of peaks
    %TODO: to modifie this to adapt for looking for just 1 peak or more
    %than 2 peaks
    %initialize onset
    p1Ind = find(onsetPoint<onsetTime_ref_peak1);
    p2Ind = find(onsetPoint<onsetTime_ref_peak2 & onsetPoint>onsetTime_ref_peak1);
    if length(p1Ind)==1
        Peak1_OnsetLat(i) = onsetPoint(p1Ind)*dx; %real onset time in msec
        Peak1_OnsetAmp(i) = onsetAmp(p1Ind); % onset amplitude in real traces (not subtracting from baseline)
        Peak1_Lat(i) = peakPoint(p1Ind)*dx; %real peak lat in msec;
        Peak1_Amp(i) = realPeakAmp(p1Ind); %real amplitude subtracting from prestim baseline
    elseif length(p1Ind)>1 %if more than one get the first one 
        Peak1_OnsetLat(i) = onsetPoint(p1Ind(1))*dx; %real onset time in msec
        Peak1_OnsetAmp(i) = onsetAmp(p1Ind(1)); % onset amplitude in real traces (not subtracting from baseline)
        Peak1_Lat(i) = peakPoint(p1Ind(1))*dx; %real peak lat in msec;
        Peak1_Amp(i) = realPeakAmp(p1Ind(1)); %real amplitude subtracting from prestim baseline
    end
    
    if length(p2Ind)==1
        Peak2_OnsetLat(i) = onsetPoint(p2Ind)*dx; %real onset time in msec
        Peak2_OnsetAmp(i) = onsetAmp(p2Ind);
        Peak2_Lat(i) = peakPoint(p2Ind)*dx; %real peak lat in msec;
        Peak2_Amp(i) = realPeakAmp(p2Ind);
   elseif length(p2Ind)>1 %if more than one get the first one 
        Peak1_OnsetLat(i) = onsetPoint(p2Ind(1))*dx; %real onset time in msec
        Peak1_OnsetAmp(i) = onsetAmp(p2Ind(1)); % onset amplitude in real traces (not subtracting from baseline)
        Peak1_Lat(i) = peakPoint(p2Ind(1))*dx; %real peak lat in msec;
        Peak1_Amp(i) = realPeakAmp(p2Ind(1)); %real amplitude subtracting from prestim baseline
    end
    
    %% plot
    nexttile
    plot(time,data); hold on;
    % peaks
    scatter(Peak1_Lat(i) +stimStart,Peak1_Amp(i)+preStimBaseline(i),'ko');hold on;
    scatter(Peak2_Lat(i) +stimStart,Peak2_Amp(i)+preStimBaseline(i),'ro');hold on;
    xlim([stimStart-20 stimStart+100])
    title(strcat('trace #', num2str(i)))
    end
end 
%% store data
PreStimBaseline = {preStimBaseline'};
Rin = {rin'};
Rs = {rs'};
Cm = {cm'};
DetectAP = {detectAP'};
NumTraces = {traceNum'};

% Peak1 table
peak1.lat = Peak1_Lat';
peak1.amp = Peak1_Amp';
peak1.onsetLat = Peak1_OnsetLat';
peak1.onsetAmp = Peak1_OnsetAmp';
Peak1 = {struct2table(peak1)};

% Peak2 table 
peak2.lat = Peak1_Lat';
peak2.amp = Peak1_Amp';
peak2.onsetLat = Peak1_OnsetLat';
peak2.onsetAmp = Peak1_OnsetAmp';
Peak2 = {struct2table(peak2)};

% All peak cell array 
AllPeakLat = {all_PeakLat}; %cell array to store all the cells 
AllPeakAmp = {all_PeakAmp};

% max slope
MaxLeftSlope = {maxLeftSlope};

% max peak amp from baseline 
MaxPeakAmp = {maxPeakAmp};

%  AUC (time integral of voltage) 
AUC_mVsec = {PSP_auc};

% definen the extra variables
Var1 = {extraVars.var1}; var1Name = extraVars.var1Name;
Var2 = extraVars.var2; var2Name = extraVars.var2Name;
Var3 = extraVars.var3; var3Name = extraVars.var3Name;


% write all the variables into a table
resTable = table(protNum,polarity,Var1,Var2,Var3,...
    PreStimBaseline,Rin,Rs,Cm,NumTraces,DetectAP,...
    MaxPeakAmp,MaxLeftSlope,AUC_mVsec,Peak1,Peak2,AllPeakLat,AllPeakAmp);
resTable.Traces = {traces'};
resTable.stimTraces = {stimTraces};
resTable.Setup = {Setup};

% write the name of the extra vars 
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

%% save 
protPath = findProtPath(cellID,protNum,meta.EphysDate,meta.base_dataPath);
save([protPath resTableName '.mat'],'resTable');

end
