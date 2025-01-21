function [resultTableDraft,traces,stimTraces,plottingStruct] = detectMinPeak(cellID,protNum,meta,polarity,monoSetup)
% Function to analyze single peak (ex. minimal stim protocal as long as
% we don't look at kinetics)
% This script will go through each trace, and detect if there are any
% single or double peaks, and will output a resultsTableDraft, which can be
% changed and saved later.
% Inputs:
% 1. cellID - use "cell1" format (no space)
% 2. protNum - use number format
% 3. meta - setting structure that has all the info about the recording day/cell
% 4. polarity - value indicating if IPSCs or EPSCs
% 5. monoSetup - structure includes predefined changable variables
% Outputs:
% 1. resultTableDraft: output a table for to be later edited per protocal.
% 2. traces: the filtered traces to be later accessed and saved

%Created by Anne - March 2019; edited by Lucas & Christine Nov-Dec 2020;
%Christine turned into function in Jan 2021
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
aoInd = find(contains(AOfields,'estim','IgnoreCase',true));
StimParams = allStimTable.(AOfields{aoInd});


% check if all the StimParams are the same for each acq, if they are, just
% store them as one single row in the stim table
if size(unique(StimParams),1) == 1
    StimParams = unique(StimParams);
end

%Parameters for the stimulation protocol (in ms):
stimStart = StimParams.delay;
searchStart = stimStart + monoSetup.searchAfterStim; %when to stop searching for peak after stim
searchStop = searchStart + monoSetup.searchWindow; %when to stop searching for peak after stim
preStim = stimStart - monoSetup.timePreStim; %time before artifact to average baseline

%find RC check params directly from the stim table
if isempty(allStimTable.RC)
    RCcheck = 0;
else
    RCcheck =1;
    pulseStart = allStimTable.RC.delay(1);
    pulseStop = pulseStart + allStimTable.RC.pulseWidth(1);
end


%% Define experiment-specific parameters

%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
dx = monoSetup.dx;

%define baseline before Rs step:
baselineStart = 1; %start computing baseline
baselineEnd = pulseStart -1; %finish computing baseline (1 ms before RC check)

%Criteria to detect events - could be modified
peakAmplitudeRef = monoSetup.peakAmplitudeRef;
doublePeakDist = monoSetup.doublePeakDist;

lowEst = monoSetup.lowEst;
highEst = monoSetup.highEst;

%% Initialize variables:
traceNum = size(traces,2);

baseline = NaN(1,traceNum);
preStimBaseline = NaN(1,traceNum);

realPeakTimeSec = NaN(1,traceNum);
finalSecPeakTime = NaN(1,traceNum);
finalSecPeakAmp = NaN(1,traceNum);
peakAmpBaselineSecPeak = NaN(1,traceNum);
peakTimefromStimSecPeak = NaN(1,traceNum);
plot_onset_point = NaN(1,traceNum); %CJL added this on 1/6/2021 -- for some reason when not initialized some of traces won't work


realPeakTimeFirst = NaN(1,traceNum);
finalFirstPeakTime = NaN(1,traceNum);
finalFirstPeakAmp = NaN(1,traceNum);
peakAmpBaselineFirstPeak = NaN(1,traceNum);
peakTimefromStimFirstPeak = NaN(1,traceNum);

twenty_percent_amp = NaN(1,traceNum);
twenty_percent_time = NaN(1,traceNum); %not really necessary to initialize them!----
seventy_percent_amp = NaN(1,traceNum);
seventy_percent_time = NaN(1,traceNum);

lin_FR_plot = cell(100,traceNum);
lin_FR_time_plot = cell(100,traceNum);
lin_FB_plot = cell(100,traceNum);
lin_FB_time_plot = cell(100,traceNum);

%Turn time into points using sampling rate dx:
%datasearchAfterStim =  searchAfterStim/dx; %2021.01.04_LV: not in use
datastimStart = stimStart/dx;
datasearchStart = searchStart/dx; %505 ms = 5050 data points
datasearchStop = searchStop/dx; %525 ms = 5250 data points
datapulseStart = pulseStart/dx;
datapulseStop = pulseStop/dx;
databaselineStart = baselineStart/dx;
databaselineEnd = baselineEnd/dx;
datapulseLength = (pulseStop-pulseStart)/dx;
datapreStim = preStim/dx;
%datasearchOnsetStop = datastimStart+(searchOnsetStop/dx);

sampling = 1:size(traces,1); %sampling points
time = sampling*dx; %in ms
filteredTraces = NaN(size(traces));
% LV_11/10/2022: now "polarity" is passed as an input and defined in master_ephys
% % Criteria for analyzing EPSCs or IPSCs:
% if contains(protName,'patch at -70mV') || contains(protName, 'min stim')
%     %mean(mean(traces(datasearchStart:datasearchStop+200,1:end)))-mean(mean(traces(datapreStim:datastimStart,1:end)))<0
%     polarity = -1; %EPSC
% elseif contains(protName,'patch at 0mV')
%     polarity = 1; %IPSC
% end
%% Find peaks:

%Define peaks trace by trace:
for i=1:traceNum %size(traces,2) %for every column in traces, that is, for every trace
    %data = traces(:,i); %original data loading
    %data = smooth(traces(:,i)); %works ok but with filter is better
    data = smooth(filter(Hd1,traces(:,i))); %smooth and filter each trace; Hd1 is custom made filter, loaded at the beggining
    
    %Compute Rs and Rm for each trace
    if RCcheck
        [rin(i), rs(i), cm(i), error(i)] = calcRs_VC(data, dx, allStimTable.RC(i,:));
    else
        rin(i) = 0;
        rs(i) = 0;
        cm(i) = 0;
        error(i) = 0;        
    end
    
    %compute baselines for each trace:
    baseline(i) = mean(data(pulseStop/dx+500:stimStart/dx-10)); %baseline between Rm step and just before stim start; more reliable than pre-pulse
    preStimBaseline(i)= mean(data(datapreStim:datastimStart)); %baseline 50 ms before stimulus
    
    %write filtered traces
    thisfilter = data;
    thisfilter(thisfilter>preStimBaseline(i)+20) = NaN;
    thisfilter(thisfilter<preStimBaseline(i)-50) = NaN;
    thisfilter(1:datastimStart-50)= NaN;
    thisfilter(datastimStart+500:end)= NaN;
    filteredTraces(:,i) = thisfilter;
    
    %Compute peak of synaptic current, EPSCs only!
    %To use this for IPSCs we should add the if polarity = +1/-1 condition:
    if polarity == -1 %if EPSC
        data = data.*-1; %need to invert traces in order to findpeaks() to work for EPSCs (positive values only)-->we'll treat the data as max, but it will be min
    end
    
    [datapeakAmplitude, datapeakTime] = findpeaks(data(datasearchStart:datasearchStop),'MinPeakDistance',10,'MinPeakProminence',1);
    peaksVector = [datapeakAmplitude, datapeakTime]; %vector with the peaks found
    
    if  isempty(peaksVector)
        detect_mono(i) = 0;
        detect_disy(i) = 0;
    else
        
        peaksVectorAmp = peaksVector(:,1);
        maxThree = maxk(peaksVectorAmp,3);
        maxPosTime = find(peaksVector==maxThree(1));
        realMaxAmp = maxThree(1);
        timeMax = peaksVector(maxPosTime,2);
        
        secMaxAmp = NaN;
        timeSecMax = NaN;
        timeThirdMax = NaN;
        if size(maxThree,1) == 2
            secMaxPosTime = find(peaksVector==maxThree(2));
            timeSecMax = peaksVector(secMaxPosTime,2);
            if abs(timeMax-timeSecMax) <= doublePeakDist
                secMaxAmp = maxThree(2);
            else
                timeSecMax = NaN;
            end
        elseif size(maxThree,1) == 3
            secMaxPosTime = find(peaksVector==maxThree(2));
            timeSecMax = peaksVector(secMaxPosTime,2);
            thirdMaxPosTime = find(peaksVector==maxThree(3));
            timeThirdMax = peaksVector(thirdMaxPosTime,2);
            if abs(timeMax-timeSecMax) <= doublePeakDist
                secMaxAmp = maxThree(2);
            elseif abs(timeMax-timeThirdMax) <= doublePeakDist
                secMaxAmp = maxThree(3);
                timeSecMax = timeThirdMax;
            end
        end
        
        if polarity == -1 %if EPSC
            data = data*-1; %data uninverted
        end
        
        %Comparison of the possible peaks with the reference amplitude:
        %(SO FAR THIS IS MY ONLY CRITERION)
        if abs(abs(realMaxAmp)-abs(preStimBaseline(i))) > peakAmplitudeRef
            firstPeakAmp = realMaxAmp;
            firstPeakTime = timeMax;
            detect_mono(i) = 1;
            if isnan(secMaxAmp) == 0 && abs(secMaxAmp) > peakAmplitudeRef+abs(preStimBaseline(i))
                if timeSecMax < timeMax
                    firstPeakAmp = secMaxAmp;
                    firstPeakTime = timeSecMax;
                    secPeakAmp = realMaxAmp;
                    secPeakTime = timeMax;
                    detect_disy(i) = 1;
                else
                    firstPeakAmp = realMaxAmp;
                    firstPeakTime = timeMax;
                    secPeakAmp = secMaxAmp;
                    secPeakTime = timeSecMax;
                    detect_disy(i) = 1;
                end
            else
                firstPeakAmp = realMaxAmp;
                firstPeakTime = timeMax;
                secPeakAmp = NaN;
                secPeakTime = NaN;
                detect_disy(i) = 0;
            end
        else
            firstPeakAmp = NaN;
            firstPeakTime = NaN;
            secPeakAmp = NaN;
            secPeakTime = NaN;
            detect_mono(i) = 0;
            detect_disy(i) = 0;
        end
        
        if polarity == -1
            %Values for firstPeakAmp and secPeakAmp (if not NaN) are converted to negative:
            firstPeakAmp = firstPeakAmp*-1;
            secPeakAmp = secPeakAmp*-1;
        end
        
        realPeakTimeSec(i) = datasearchStart + secPeakTime-1;
        finalSecPeakTime(i) = realPeakTimeSec(i)*dx;
        finalSecPeakAmp(i) = secPeakAmp;
        peakAmpBaselineSecPeak(i) = finalSecPeakAmp(i)-preStimBaseline(i);
        peakTimefromStimSecPeak(i) = finalSecPeakTime(i) - stimStart;
        
        realPeakTimeFirst(i) = datasearchStart + firstPeakTime-1; %I need to store it for kinetics
        finalFirstPeakTime(i) = realPeakTimeFirst(i)*dx;
        finalFirstPeakAmp(i) = firstPeakAmp;
        peakAmpBaselineFirstPeak(i) = finalFirstPeakAmp(i)-preStimBaseline(i);
        peakTimefromStimFirstPeak(i) = finalFirstPeakTime(i) - stimStart;
        
        %% onset
        %Now we are going to find the onset of the 1st peak (or real peak):
        %We calculate the 20% amplitude of the peak and make it relative to the
        %data context; then we do a walk-back (OP1) or an interpolation (OP2):
        if detect_mono(i) == 1
            %------------------option 1: walk back/walk forward to find 20% and 70% of peak amp in real data------------------
            %    to find 20% value in data: walk back from peak time until the stim (dataSearchStart)
            ideal_twenty_percent = lowEst*(data(realPeakTimeFirst(i))-data(datasearchStart)); %real 20% value
            twenty_percent_amp(i) = ideal_twenty_percent+data(datasearchStart); %the amplitude of that ideal 20%
            k = realPeakTimeFirst(i);
            k2=k-length(datasearchStart:realPeakTimeFirst(i))+1; %k-length(data_amp);
            foundk = 0;
            
            while k > k2
                if abs(twenty_percent_amp(i)-data(k)) > 1 %look for values in the data that are less than 5 pA away from ideal 20%; 5 pA can be modified
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
            ideal_seventy_percent = highEst*(data(realPeakTimeFirst(i))-data(datasearchStart)); %real 20% value
            seventy_percent_amp(i) = ideal_seventy_percent+data(datasearchStart); %the amplitude of that ideal 20%
            j = k;
            j2 = j+length(twenty_percent_time(i)/dx:realPeakTimeFirst(i))-1;
            
            foundj = 0;
            while j < j2 
                if abs(seventy_percent_amp(i)-data(j)) > 1 %look for values in the data that are less than 5 pA away from ideal 70%; 5 pA can be modified
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
            %----------------------------------------OPTION 2-----------------------------------------
%             
%             % now we interpolate data (only between datasearchStart and the first real peak
%             r = 10000; %(how many times more points do we want to add to the sampling data
%             x = datasearchStart:realPeakTimeFirst(i); % x axis data points (sample)
%             v = data([datasearchStart:realPeakTimeFirst(i)]); % the cooresponding values of the sample data points
%             xq = datasearchStart:1/r:realPeakTimeFirst(i); % increase the sampling rate by r
%             interp_data = interp1(x,v,xq); % this is the interpolation (default style is linear, should be ok here)
%             
%             % ---- to find 20% ----
%             ideal_twenty_percent = lowEst*(data(realPeakTimeFirst(i))-preStimBaseline(i)); %real 20% value
%             twenty_percent_amp(i) = ideal_twenty_percent+preStimBaseline(i); %the amplitude of that ideal 20%
%             
%             % the "real" twenty_percent_amp in the interpolated data
%             logical_sample = abs(twenty_percent_amp(i)-interp_data) < 3;
%             onlyT_data = interp_data(logical_sample); % only the data that is true for the previous logical statement
%             twenty_percent_amp(i) = onlyT_data(length(onlyT_data)); %the last of the logical True answers is the equivalent to the walk back)
%             
%             % to find the twenty_percent_time(i)
%             z = (find(interp_data==twenty_percent_amp(i))); % in case there are more than one data point has the exact value of the twenty_per_cent amplitudate
%             z = z(length(z)); %  z is now the last index that matches the amplitude of the twenty_percent_amp in the interpolated data
%             k = round(datasearchStart+(z/r)); % for later use
%             twenty_percent_time(i) = (datasearchStart+(z/r))*dx; % now convert to time
%             
%             
%             % ---- to find 70% ----
%             ideal_seventy_percent = highEst*(data(realPeakTimeFirst(i))-preStimBaseline(i)); %real 70% value
%             seventy_percent_amp(i) = ideal_seventy_percent+preStimBaseline(i); %the amplitude of that ideal 70%
%             
%             logical_sample_70 = abs(seventy_percent_amp(i)-interp_data) < 2;
%             onlyT_data_70 = interp_data(logical_sample_70); % only the data that is true for the previous logical statement
%             seventy_percent_amp(i) = onlyT_data_70(length(onlyT_data_70)); %the first of the logical True answers is the equivalent to the walk forward)
%             
%             y = (find(interp_data==seventy_percent_amp(i))); % in case there are more than one data point has the exact value of the twenty_per_cent amplitudate
%             y = y(1); %  y is now the first index that matches the amplitude of the seventy_percent_amp in the interpolated data
%             j = round(datasearchStart+(y/r)); % for later use
%             seventy_percent_time(i) = j*dx; % now convert to time
            %
            %
            %
            %         %---------------------------------------END of OPTION 2---------------------------------------
            %         x1 = k:round(realPeakTimeFirst(i)); %or finalFirstPeakTime(i)/dx
            x1 = k:j;
            if length(x1)<=1
                onsetTime(i) = searchStart-stimStart + 1; %default is 1 second after stim if no onset could be calculated
                onsetAmp(i) = data((searchStart+1)./dx);
                twenty_percent_amp(i) = NaN;
                twenty_percent_time(i) = NaN;
            else
                FR = fit(x1', data(x1), 'poly1'); %same as using polyfit(x1',data(x1),1) --> linear fit
                CR = coeffvalues(FR);
                lin_FR_vec = k:j;
                lin_FR_dist = 100-length(lin_FR_vec); %I want 100 pts lines, so this is the number of points to add
                lin_FR = CR(1)*(k-lin_FR_dist:j)+CR(2); %this has 100 pts
                lin_FR_plot{i} = lin_FR;
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
                
                %NEW2 OK:
                fit_onset_amp = lin_FR(intersect_time); %the closest value in lin_FR to the real data
                closest_to_intersect = NaN(1,length(k-lin_FR_dist:realPeakTimeFirst(i)));
                %%empty vector to fill w/ closest values in real data to the amplitude of
                %%intersect of lin_FR and lin_FB
                %
                %t_vector = k-89:k;
                %time points to look for those values
                t_vector = lin_FR_time(intersect_time)-(abs(k-lin_FR_time(intersect_time))):k; %the time in points of the
                %intercept
                if length(t_vector) < 3
                    onsetAmp(i) = data(k);
                    onsetTime(i) =  k*dx-stimStart;
                else
                    for j = 1:length(t_vector)
                        closest_to_intersect(j) = min(abs(fit_onset_amp-data(t_vector(j))));
                    end
                    
                    %among values in closest_to_intersect I want theclosest to 0 with highest index:
                    [min_val,min_indices] = mink(closest_to_intersect,3);
                    onset_point(i) = t_vector(max(min_indices)); %point time in which the onset happens
                    if onset_point(i) < (stimStart + monoSetup.searchAfterStim + 1)/dx
                        onset_point(i) =(stimStart + monoSetup.searchAfterStim + 1)/dx;
                    end
                    plot_onset_point(i) = onset_point(i)*dx; %converted to sweep time (ms)
                    onsetAmp(i) = data(onset_point(i));
                    onsetTime(i)= plot_onset_point(i)-stimStart; %real onset latency
                end
            end
            
        else
            onsetTime(i) = NaN;
            onsetAmp(i) = NaN;
            twenty_percent_amp(i) = NaN;
            twenty_percent_time(i) = NaN;
        end
    end
    
end

%% Store results:

results_array = [preStimBaseline;rin;rs;cm;detect_mono;onsetTime;onsetAmp;finalFirstPeakTime;finalFirstPeakAmp;...
    detect_disy;finalSecPeakTime;finalSecPeakAmp];
results_array = results_array';
resultTableDraft = array2table(results_array,'VariableNames',...
    {'preStimBaseline','rin','rs','cm','detect_mono','onsetTime','onsetAmp','finalFirstPeakTime','finalFirstPeakAmp',...
    'detect_disy','finalSecPeakTime','finalSecPeakAmp'});

% Required for the app:
plottingStruct.filteredTraces = filteredTraces;
plottingStruct.stimStart = stimStart;
% plottingStruct.onsetAmp = onsetAmp;
%% Plot responses:
%All the traces together (will do this inside the loop w/ subplot, but for quick

%computer filtered traces
% filteredTraces = NaN(size(traces));
% for i=1:traceNum
%     filteredTraces(:,i) = smooth(filter(Hd1,traces(:,i)));
% end
% filteredTraces = filteredTraces';

% %%
% for i=1:traceNum
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     plot(time, filteredTraces(i,:)); %plot filtered response
%     xlabel('Time (ms)');
%     ylabel('Amplitude (pA)');
%     figure_label = ['Trace ' num2str(i)];
%     title(figure_label);
%     if detect_mono(i) == 1 && detect_disy(i) == 1 %mono and disynaptic
%         hold on; plot(finalFirstPeakTime(i),finalFirstPeakAmp(i),'bo');
%         hold on; plot(finalSecPeakTime(i),finalSecPeakAmp(i),'ro');
%         hold on; plot(plot_onset_point(i),onsetAmp(i),'c^');
%         hold on; plot(twenty_percent_time(i),twenty_percent_amp(i),'ko');
%         hold on; plot(lin_FR_time_plot{i}*dx,lin_FR_plot{i},'c-');
%         hold on; plot(lin_FB_time_plot{i}*dx,lin_FB_plot{i},'r-');
%         xlim([stimStart-20 stimStart+150]);
%         ylim([preStimBaseline(i)-75 preStimBaseline(i)+20]);
%         ax = gca;
%         ax.XColor = 'yellow';
%     elseif detect_mono(i) == 1 && detect_disy(i) == 0 %monosynaptic
%         hold on; plot(finalFirstPeakTime(i),finalFirstPeakAmp(i),'bo');
%         hold on; plot(plot_onset_point(i),onsetAmp(i),'c^');
%         hold on; plot(twenty_percent_time(i),twenty_percent_amp(i),'ko');
%         hold on; plot(lin_FR_time_plot{i}*dx,lin_FR_plot{i},'c-');
%         hold on; plot(lin_FB_time_plot{i}*dx,lin_FB_plot{i},'r-');
%         xlim([stimStart-20 stimStart+150]);
%         ylim([preStimBaseline(i)-75 preStimBaseline(i)+20]);
%         ax = gca;
%         ax.XColor = 'green';
%
%     else
%         xlim([stimStart-20 stimStart+150]);
%         ylim([preStimBaseline(i)-75 preStimBaseline(i)+20]);
%         ax = gca;
%         ax.XColor = 'red';
%     end
% end

end
