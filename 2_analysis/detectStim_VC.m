function [resultTableDraft,traces,stimTraces] = detectStim_VC(cellID,protNum,meta,polarity,Setup,stimType)

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
%Updated to simplify loading data and stim & updated calcRs_VC integration by Christine on 3/10/23

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

%Criteria to detect events - could be modified 
% stimulation protocols:
peakAmplitudeRef = Setup.peakAmplitudeRef;
onsetTime_ref_EPSC = Setup.onsetTime_ref_EPSC;
onsetTime_ref_IPSC = Setup.onsetTime_ref_IPSC;
fallTimefromPeak_ref = Setup.fallTimefromPeak_ref;
totalCharge_ref_IPSC = Setup.totalCharge_ref_IPSC;
totalCharge_ref_EPSC = Setup.totalCharge_ref_EPSC;
fallTimeMax = Setup.fallTimeMax; %max possible return to baseline
% for estimating onset 
lowEst = Setup.lowEst;
highEst = Setup.highEst;

%% Initialize variables:
traceNum = size(traces,2);

baseline = NaN(1,traceNum);
preStimBaseline = NaN(1,traceNum);

realpeakTime = NaN(1,traceNum);
peakTime = NaN(1,traceNum);
peakAmplitude = NaN(1,traceNum);
peakAmplitudeBaseline = NaN(1,traceNum);
peakTimefromStim = NaN(1,traceNum);

twenty_percent_amp = NaN(1,traceNum);
twenty_percent_time = NaN(1,traceNum);

seventy_percent_amp = NaN(1,traceNum);
seventy_percent_time = NaN(1,traceNum);

lin_FR_plot = cell(100,traceNum);
lin_FR_time_plot = cell(100,traceNum);
lin_FB_plot = cell(100,traceNum);
lin_FB_time_plot = cell(100,traceNum);

%initialize all storage variable as NaN:
detect = zeros(1,traceNum); 
riseTime = NaN(1,traceNum); 
onsetTime = NaN(1,traceNum);
fallTimefromPeak = NaN(1,traceNum);
tauFall = NaN(1,traceNum);
rsquareMonoFit = NaN(1,traceNum);
tauFallDoubleFirst = NaN(1,traceNum);
tauFallDoubleSecond = NaN(1,traceNum);
rsquareDoubleFit = NaN(1,traceNum);


%Turn time into points using sampling rate dx:
datastimStart = stimStart/dx; 
datasearchStart = searchStart/dx; 
datasearchStop = searchStop/dx; 
datapreStim = preStim/dx;

%----

sampling = 1:size(traces,1); %sampling points
time = sampling*dx;

%% Analyze traces:

figure('units','normalized','outerposition',[0 0 1 1]); %creates fig that will be filled with subplot

%Trace by trace analysis:
for i=1:size(traces,2) %specifies that we refer to 2nd col: 50 (# of traces) and not to the 1st (number of time points)
    %data = traces(:,i);
    data = smooth(filter(Hd1,traces(:,i))); %smooth and filter each trace; Hd1 is custom made filter, loaded at the beggining
    
    %compute baselines for each trace
     
    baseline(i) = mean(data(1:pulseStart)); % baseline before Rm step
    preStimBaseline(i)= mean(data(datapreStim:datastimStart)); %baseline 50 ms before stimulus
   
    
    %compute peak of synaptic current, taking into account if EPSC or IPSC:
    if polarity == 1
        [datapeakAmplitude, datapeakTime] = max(data(datasearchStart:datasearchStop));
        %IPSCs: no need to pick more than one peak
    else
        [datapeakAmplitude, datapeakTime] = min(data(datasearchStart:datasearchStop));
        %EPSCs: find all the peaks (not just max) and choose the one
        %closest to the stim
    end
    
    realpeakTime(i) = datasearchStart + datapeakTime-1;
    peakTime(i) = realpeakTime(i)*dx;
    peakAmplitude(i) = datapeakAmplitude;
    peakAmplitudeBaseline(i) = peakAmplitude(i)-preStimBaseline(i);
    peakTimefromStim(i) = peakTime(i) - stimStart;
    
    
    %Compute Rs and Rm for each trace
    if RCcheck
        [rin(i), rs(i), cm(i), error(i)] = calcRs_VC(traces(:,i), dx, allStimTable.RC(i,:));
    else 
        rin(i) = 0;conda
         rs(i) = 0;
         cm(i) = 0;
       error(i) = 0;
         
    end 
    
    %Compute onset time of synaptic current
    
        %------------------option 1: walk back/walk forward to find 20% and 70% of peak amp in real data------------------
        %    to find 20% value in data: walk back from peak time until the stim (dataSearchStart)
            ideal_twenty_percent = 0.2*(data(realpeakTime(i))-data(datasearchStart)); %real 20% value
            twenty_percent_amp(i) = ideal_twenty_percent+data(datasearchStart); %the amplitude of that ideal 20%
            k = realpeakTime(i);
            k2=k-length(datasearchStart:realpeakTime(i))+1; %k-length(data_amp);
            foundk = 0;
    
            while k > k2
                if abs(twenty_percent_amp(i)-data(k)) > 8 %look for values in the data that are less than 5 pA away from ideal 20%; 5 pA can be modified
                    k = k-1;
                    foundk = 0;
                else
                    foundk = 1;
                    %k=k;
                    break
                end
                
            end
            twenty_percent_amp(i) = data(k);
            twenty_percent_time(i) = k*dx;
            
            
            %to find 70% value in data: walk forward from the 20% time until the peak time
            ideal_seventy_percent = 0.7*(data(realpeakTime(i))-data(datasearchStart)); %real 20% value
            seventy_percent_amp(i) = ideal_seventy_percent+data(datasearchStart); %the amplitude of that ideal 20%
            j = k;%twenty_percent_time(i)/dx;
            j2 = j+length(twenty_percent_time(i)/dx:realpeakTime(i))-1;
            foundj = 0;
            while j < j2 %j <= j2
                if abs(seventy_percent_amp(i)-data(j)) > 8 %look for values in the data that are less than 5 pA away from ideal 20%; 5 pA can be modified
                    j = j+1;
                    foundj = 0;
                else
                    foundj = 1;
                    %j=j;
                    break
                end
            end
    
            seventy_percent_amp(i) = data(j);
            seventy_percent_time(i) = j*dx;
           
            
       %---------------------------------------end of option 1 ------------------------------------------
    
    
    
%     % ---------------------------------------- option 2  ---------------------------------------
%     % set up interpolation between datasearch onset and peak
%     r = 500; %(how many times more points do we want to add to the sampling data
%     x = datasearchStart:realpeakTime(i); % x axis data points (sample)
%     v = data([datasearchStart:realpeakTime(i)]); % the cooresponding values of the sample data points
%     xq = datasearchStart:1/r:realpeakTime(i); % increase the sampling rate by r
%     interp_data = interp1(x,v,xq); % this is the interpolation (default style is linear, should be ok here)
%     
%     % ---- to find 20% ----
%     % define ideal 20% amp to look for
%     ideal_twenty_percent = 0.2*(data(realpeakTime(i))-data(datasearchStart)); %real 20% value
%     twenty_percent_amp(i) = ideal_twenty_percent+data(datasearchStart); %the amplitude of that ideal 20%
%     
%     % the "real" twenty_percent_amp in the interpolated data
%     logical_sample = abs(twenty_percent_amp(i)-interp_data) < 0.5;
%     onlyT_data_20 = interp_data(logical_sample); % only the data that is true for the previous logical statement
%     twenty_percent_amp(i) = onlyT_data_20(length(onlyT_data_20)); %the last of the logical True answers is the equivalent to the walk back)
%     
%     % to find the twenty_percent_time(i)
%     z = (find(interp_data==twenty_percent_amp(i))); % in case there are more than one data point has the exact value of the twenty_per_cent amplitudate
%     z = z(length(z)); %  z is now the last index that matches the amplitude of the twenty_percent_amp in the interpolated data
%     k = round(datasearchStart+(z/r)); % for later use -- conversion to k is not perfect
%     %twenty_percent_time(i) = (datasearchStart+(z/r))*dx; % now convert to time
%     twenty_percent_time(i) = k*dx;
%     
%     % ---- to find 70% ----
%     
%     ideal_seventy_percent = 0.7*(data(realpeakTime(i))-data(datasearchStart)); %real 60% value
%     seventy_percent_amp(i) = ideal_seventy_percent+data(datasearchStart); %the amplitude of that ideal 20%
%     
%     logical_sample_70 = abs(seventy_percent_amp(i)-interp_data) < 0.5;
%     onlyT_data_70 = interp_data(logical_sample_70); % only the data that is true for the previous logical statement
%     seventy_percent_amp(i) = onlyT_data_70(1); %the first of the logical True answers is the equivalent to the walk forward)
%     
%     y = (find(interp_data==seventy_percent_amp(i))); % in case there are more than one data point has the exact value of the twenty_per_cent amplitudate
%     y = y(1); %  y is now the first index that matches the amplitude of the seventy_percent_amp in the interpolated data
%     j = round(datasearchStart+(y/r)); % for later use
%     seventy_percent_time(i) = j*dx; % now convert to time
%     % ---------------------------------------- END of option 2  ---------------------------------------
    
    
    x1 = k:round(realpeakTime(i)); %or finalFirstPeakTime(i)/dx
%     x1 = k:j;
     if length(x1) >= 2 %requirement of 2 points or more to fit rise; if id doesn't find points everything crashes (this happens when responses are abolished)
        
        FR = fit(x1', data(x1), 'poly1'); %same as using polyfit(x1',data(x1),1) --> linear fit
        CR = coeffvalues(FR);
        %lin_FR_vec = k:realpeakTime(i);
        lin_FR_vec = k:j;
        lin_FR_dist = 100-length(lin_FR_vec); %I want 100 pts lines, so this is the number of points to add
        %lin_FR = CR(1)*(k-lin_FR_dist:realpeakTime(i))+CR(2); %this has 100 pts
        lin_FR = CR(1)*(k-lin_FR_dist:j)+CR(2); %this has 100 pts
        lin_FR_plot{i} = lin_FR;
        %lin_FR_time = realpeakTime(i)-length(lin_FR)+1:realpeakTime(i); %and this is the time axis
        lin_FR_time = j-length(lin_FR)+1:j; %and this is the time axis
        lin_FR_time_plot{i} = lin_FR_time;
        x2 = pulseStop/dx+500:stimStart/dx-10;
        FB = fit(x2', data(x2),'poly1');
        CB = coeffvalues(FB);
        lin_FB = CB(1)*(datasearchStart-5:datasearchStart+length(lin_FR)-6)+CB(2);
        lin_FB_plot{i} = lin_FB;
        lin_FB_time = datasearchStart-5:datasearchStart+length(lin_FB)-6;
        lin_FB_time_plot{i} = lin_FB_time;
        
        [intersect_val, intersect_time] = min(abs(lin_FR-lin_FB));
        
        fit_onset_amp = lin_FR(intersect_time); %the closest value in lin_FR to the real data
        %closest_to_intersect = NaN(1,length(k-lin_FR_dist:realpeakTime(i)));
        closest_to_intersect = NaN(1,length(k-lin_FR_dist:j));
        %%empty vector to fill w/ closest values in real data to the amplitude of
        %%intersect of lin_FR and lin_FB
        %
        %t_vector = k-89:k;
        %time points to look for those values
        
        if lin_FR_time(intersect_time)-(abs(k-lin_FR_time(intersect_time))) > datasearchStart % to make sure t_vector starts no earlier than datasearchStart
            t_vector = lin_FR_time(intersect_time)-(abs(k-lin_FR_time(intersect_time))):k; %the time in points of the intercept
        else
            t_vector = datasearchStart : k;
        end
        %intercept
        if length(t_vector) < 3
            onsetAmp(i) = data(k);
            onsetTime(i) =  k*dx-stimStart;%added '-stimStart' (LV_01132021)
            plot_onset_point(i) = k*dx; %added (LV_01132021)
            onset_point(i)=plot_onset_point(i)/dx; %added (LV_01132021)
            dataonsetTime = round(onset_point(i));
            
            
        else
            
            for k = 1:length(t_vector)
                closest_to_intersect(k) = min(abs(fit_onset_amp-data(t_vector(k))));
            end
         
            %among values in closest_to_intersect I want theclosest to 0 with highest index:
            [min_val,min_indices] = mink(closest_to_intersect,3);
            onset_point(i) = t_vector(max(min_indices)); %point time in which the onset happens
            
            plot_onset_point(i) = onset_point(i)*dx; %converted to sweep time (ms)
            onsetAmp(i) = data(onset_point(i));
            onsetTime(i)= plot_onset_point(i)-stimStart; %real onset latency
            dataonsetTime = round(onset_point(i));
        end
            
            
            %compute fall time to baseline (10% of amplitude)
            r =  round(peakTime(i)/dx);
            numRow = 0;
            if polarity == 1
                tenPercentAmplitude = (preStimBaseline(i)+(0.1*peakAmplitudeBaseline(i)));
                while (numRow < 3 && r<size(traces,1))
                    if data(r) < tenPercentAmplitude
                        numRow = numRow + 1; % 3 consecutive points smaller than 20% peak
                    else
                        numRow = 0;
                    end
                    r = r+1;
                    if r>=size(traces,1)
                        disp('no synaptic current fall detected');
                    end
                end
            else
                tenPercentAmplitude = (preStimBaseline(i)+(0.1*peakAmplitudeBaseline(i)));
                while (numRow < 3 && r<size(traces,1))
                    if data(r) > tenPercentAmplitude
                        numRow = numRow + 1;
                    else
                        numRow = 0;
                    end
                    r = r+1;
                    if r>=size(traces,1)
                        disp('no synaptic current fall detected');
                    end
                end
            end
            
            fallTime(i) = (r-3)*dx;
            fallTimefromPeak(i) = fallTime(i)-peakTime(i);
            datafallTime = round(fallTime(i)/dx);
            fallAmplitude(i) = data(datafallTime);
           
                       
            %compute charge
            currentTrace = data(dataonsetTime:datafallTime)-preStimBaseline(i);
            totalCharge(i) = trapz(currentTrace)/10000;
        
        
    elseif length(x1) < 2
        plot_onset_point(i) = NaN; %converted to sweep time (ms)
        onsetAmp(i) = NaN;
        onsetTime(i)= NaN; %real onset latency
        fallTime(i) = NaN;
        fallTimefromPeak(i) = NaN;
        datafallTime = NaN;
        fallAmplitude(i) = NaN;
        currentTrace = NaN;
        totalCharge(i) = NaN;
     end
    
    
    %Criteria for defining response or not:
    %2019.12.15_LV: added norise == 0 as a condition to evaluate
    %2020.03.02_LV: changed peakTime(i)-realonsetTime(i)<20 to <50 for polarity==1, since some responses are slower at +10 mV
    if (polarity == 1 && peakAmplitudeBaseline(i)>peakAmplitudeRef &&... %for IPSC
            onsetTime(i)<onsetTime_ref_IPSC && ...
            fallTimefromPeak(i)>fallTimefromPeak_ref &&...
            abs(totalCharge(i))>abs(totalCharge_ref_IPSC) &&...
            peakTime(i)-onset_point(i)<50) || (polarity == -1 && ...
            peakAmplitudeBaseline(i)<-1*peakAmplitudeRef &&...
            onsetTime(i)<onsetTime_ref_EPSC &&...
            fallTimefromPeak(i)> fallTimefromPeak_ref &&...
            abs(totalCharge(i))> abs(totalCharge_ref_EPSC) &&...
            peakTime(i)-onset_point(i)<20)
        detect(i) = 1;
    else
        detect(i) = 0;
    end
    
    if detect(i) == 1
        %compute 10-90% rise time
        tenPercentAmplitude = (preStimBaseline(i)+(0.1*peakAmplitudeBaseline(i)));
        ninetyPercentAmplitude = (preStimBaseline(i)+(0.9*peakAmplitudeBaseline(i)));
        risingTrace = data(datasearchStart:realpeakTime(i));
        [risingTwentyPercentAmp, risingTwentyPercentTime] = min(abs(risingTrace-tenPercentAmplitude));
        [risingNinetyPercentAmp, risingNinetyPercentTime] = min(abs(risingTrace-ninetyPercentAmplitude));
        riseTime(i) = (risingNinetyPercentTime-risingTwentyPercentTime)*dx;
        
        %compute tau of falling phase
        if  fallTimefromPeak(i) < fallTimeMax % exclude abnoraml fallTime and tau due to multiple peaks
            fallingPhase = data(realpeakTime(i):datafallTime)-preStimBaseline(i);
            x = 1:length(fallingPhase);
            [fitfallingPhase, gof] = fit(x',fallingPhase,'exp1');
            realfitfallingPhase = squeeze(fitfallingPhase(x)+preStimBaseline(i));
            realfitfallingPhase_all{i,:} = realfitfallingPhase;
            datatauFall = (1/(fitfallingPhase.b*-1));
            tauFall(i) = datatauFall*dx;
            rsquareMonoFit(i) = gof.rsquare;
            [fitfallingPhaseDouble, gof] = fit(x',fallingPhase,'exp2');
            realfitfallingPhaseDouble = squeeze(fitfallingPhaseDouble(x)+preStimBaseline(i));
            realfitfallingPhaseDouble_all{i,:} = realfitfallingPhaseDouble;
            datatauFallDoubleFirst = (1/(fitfallingPhaseDouble.b*-1));
            datatauFallDoubleSecond = (1/(fitfallingPhaseDouble.d*-1));
            tauFallDoubleFirst(i) = datatauFallDoubleFirst*dx;
            tauFallDoubleSecond(i) = datatauFallDoubleSecond*dx;
            rsquareDoubleFit(i) = gof.rsquare;
        end
    end
    
    
    
    % Subplots to have a quick visualization of traces:
    if traceNum <= 5
        subpRow = 1;
        subpCol = traceNum;
    elseif (5 < traceNum) && (traceNum <= 10)
        subpRow = 2;
        subpCol = 5;
    elseif (10 < traceNum) && (traceNum <= 50)
        subpRow = 5;
        subpCol = 10;
    elseif traceNum > 50
        subpRow = 10;
        subpCol = 10;
    end
    
    %plot each trace with corresponding marks
    subplot1 = subplot(subpRow,subpCol,i); %2020.01.01_LV: subplot(5,10,i) just for 50 traces
    %figure('units','normalized','outerposition',[0 0 1 1]);
    plot(time, smooth(filter(Hd1,traces(:,i))),'Parent',subplot1); %plot response
    xlabel('Time (msec)');
    ylabel('Amplitude (pA)'); hold on;
    xlim([stimStart-20 stimStart+150]);
    figure_label = ['Trace ' num2str(i)];
    title(figure_label);
    hold on; plot(peakTime(i),peakAmplitude(i),'o','Parent',subplot1);
    %hold on; plot(fallTime(i),fallAmplitude(i),'o','Parent',subplot1);
    hold on; plot(onsetTime(i),onsetAmp(i),'x','Parent',subplot1);
    
    if detect(i) == 1
        if fallTimefromPeak(i)<fallTimeMax
            hold on; plot(peakTime(i):0.1:fallTime(i),realfitfallingPhase,'r--','LineWidth',2,'Parent',subplot1);
            hold on; plot(peakTime(i):0.1:fallTime(i),realfitfallingPhaseDouble,'g--','LineWidth',2,'Parent',subplot1);
            hold on; plot(peakTime(i),peakAmplitude(i),'o','Parent',subplot1);
            %hold on; plot(fallTime(i),fallAmplitude(i),'o','Parent',subplot1);
            hold on; plot(onsetTime(i),onsetAmp(i),'x','Parent',subplot1);
        end
    else
        ax = gca;
        ax.XColor = 'red';
    end
    %xlim(subplot1, [stimStart-20 stimStart+fallTimeMax+50]);
    xlim([stimStart-20 stimStart+150]);
    if polarity == 1
        ylim(subplot1, [preStimBaseline(i)-20 preStimBaseline(i)+peakAmplitudeBaseline(i)+40]);
    else
        ylim(subplot1, [preStimBaseline(i)+peakAmplitudeBaseline(i)-20 preStimBaseline(i)+40]);
    end
    cellType = values(meta.mapCellType,{cellID});
    sgtitle([meta.EphysDate ':' cellID '-' cellType{1} ' - ProtNum' num2str(protNum)])
end
%% Store results

results_array = [preStimBaseline; rin; rs; cm; detect; peakAmplitudeBaseline; peakTimefromStim; onsetTime;...
    totalCharge; riseTime; fallTimefromPeak; tauFall;rsquareMonoFit];%; rsquareMonoFit; tauFallDoubleFirst; tauFallDoubleSecond;...
    %rsquareDoubleFit];

results_array = results_array';

resultTableDraft = array2table(results_array,'VariableNames',...
     {'preStimBaseline','rin','rs','cm','detect','peakAmplitudeBaseline','peakTimefromStim','onsetTime',...
     'totalCharge','riseTime','fallTimefromPeak','tauFall','rsquareMonoFit'});%,'tauFallDoubleFirst',...


end