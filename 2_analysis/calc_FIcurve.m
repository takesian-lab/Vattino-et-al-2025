function [spikes_table,Rm,Rs,rheoBase] = calc_FIcurve(pathname,plotPause)
%this function calculates all the spike properties, and Rs and Rm for each
%trace in the IV curve. analysis code written by Lucas; adapted into
%function by Christine

%Will output a table with results from each trace, and Rm and Rs calcalated
%from -30pA current injection steps;

%search for TODO.
if nargin == 1
  plotPause = 1; %Defult for plot is yes;
end


%% load data traces and stim parameters 
[dataTraces,fileName] = loadData(pathname);
[stimParamsStruct,~] = loadStim(pathname);
stimTable = stimParamsStruct.intrinsic;

%% set stim parameters 
dx = 0.1;
pulseStart = stimTable.delay(1); %in ms; time of injected current
pulseLength = stimTable.pulseWidth(1); %in ms; duration of the injected pulse
pulseEnd = pulseStart + pulseLength; %in ms; end of injected current

datapulseStart = pulseStart/dx; %time of injected current in points
datapulseEnd = pulseEnd/dx; %end of injected current in points
stimAmpUsed = stimTable.amplitude;

% if this is a full IV curve (have both stim amp values smaller and greater
% than o:
if ~isempty(find(stimAmpUsed>0, 1)) && ~isempty(find(stimAmpUsed<0, 1))
    % but if there is no 0
    if isempty(find(stimAmpUsed==0, 1))
        zero_index = find(stimAmpUsed == min(stimAmpUsed(stimAmpUsed>0))); %the current first index with positive value should be the new zero
        stimAmpUsed = [stimAmpUsed(1:zero_index-1);0;stimAmpUsed(zero_index:end)];
    end
end
%% find rheoIndex

pulseIndex = [datapulseStart:datapulseEnd];
maxDuringStim = max(dataTraces(pulseIndex,:)); % the max voltage (mV) during the stim period for each trace;
indexes = find(maxDuringStim>0); % index the traces that have an AP (or max voltage is greater than 0; this will first return all the indexes that have APs
if ~isempty(indexes)
    rheoIndex = indexes(1);
    rheoBase = stimAmpUsed(rheoIndex);
else 
    rheoIndex = 1;
    rheoBase = [];
end 
%% RC check using -30pA injection
RC_check_data = dataTraces(:,stimAmpUsed==30); % = data(:,index_RC_sweep)
RC_check_Stimtble = stimTable(stimAmpUsed==30,:);
[Rm, Rs, Cm, error] = calcRs_CC(RC_check_data,dx,RC_check_Stimtble);
% [Rm, Rs, Cm, error] = deal(nan)

%% set values
datapulseStart = pulseStart/dx; %time of injected current in points
datapulseEnd = pulseEnd/dx; %end of injected current in points
%dataRCpulseStart = RCpulseStart/dx;

%[d,si]=abfload([pathname fileName]);
%data = squeeze(d);
sampling = 1:size(dataTraces,1); %vector containing with as many columns as sampling points
time = sampling*dx;


%I_inject = first_step:increase_step:last_step; %moved above to include
%non-hSGNs option
indexZero = find(stimAmpUsed == 0); %sweep where curr inj = 0


%traces_to_analyze = obs_threshold:length(I_inject); %traces that showed APs, starting from the one that's the threshold

% Rs_MOhm = nan(1,size(data,2));
% Rm_MOhm = nan(1,size(data,2));
%pulseAmp = nan(1,size(data,2));

baseline = nan(1,size(dataTraces,2));

dataSAG = nan(1,size(dataTraces,2));
dataminTime = nan(1,size(dataTraces,2));
datasteadyTime = nan(1,size(dataTraces,2));
dataminAmplitude = nan(1,size(dataTraces,2));
datasteadyAmplitude = nan(1,size(dataTraces,2));

datapeakAmplitude = cell(size(dataTraces,2),1);
datapeakTime = cell(size(dataTraces,2),1);

datathresholdTime = cell(size(dataTraces,2),1);
datathresholdAmplitude = cell(size(dataTraces,2),1);

datatroughAmplitude = cell(size(dataTraces,2),1);
datatroughTime = cell(size(dataTraces,2),1);

dataAPheight = cell(size(dataTraces,2),1);
dataISI = cell(size(dataTraces,2),1);

datahalfwidth = cell(size(dataTraces,2),1);

datahalftime1 = cell(size(dataTraces,2),1);
datahalftime2 = cell(size(dataTraces,2),1);

spike_result = nan(length(stimAmpUsed),16); %based on the number of output
%variables that I'm currently using

%% F-I curve
%Calculate spike frequency, spike thrshold and and spike number for each
%current injection trace

for i = 1:size(dataTraces,2) %obs_threshold:obs_threshold+1
    trace_data = smooth(dataTraces(:,i));
    baseline(i) = mean(trace_data(1:99)); %baseline before RC check
    %pulseAmp(i) = I_inject(i); %amplitude of pulse
    spike_result(i,1) = stimAmpUsed(i);
    spike_result(i,2) = baseline(i);
    
    if i < rheoIndex
        spike_result(i,3) = NaN;
        
        if i < indexZero %I want SAG for traces below I_inject == 0
            [minAmplitude, minTime] = min(trace_data(datapulseStart:datapulseStart+3000)); %min mV in the 1st ~third of the step
            [steadyAmplitude, steadyTime] = min(trace_data(datapulseEnd-250:datapulseEnd-1)); %min mV in the last 25 ms
            
            %if dataSAG < 0, that means that we are in the opposite
            %situation; the sign matters in this case
            dataSAG(i) = steadyAmplitude-minAmplitude; %withouth considering tau, something that I should add later (see Allen whitepaper)
            
            minTime = datapulseStart-1 + minTime;
            steadyTime = datapulseEnd-250 + steadyTime;
            dataminTime(i) = minTime;
            datasteadyTime(i) = steadyTime;
            
            dataminAmplitude(i) = trace_data(minTime);
            datasteadyAmplitude(i) = trace_data(steadyTime);
            
            
            spike_result(1,3) = dataSAG(i);
        end
        for j = 4:size(spike_result,2)
            spike_result(i,j) = NaN;
        end
    elseif i >=  rheoIndex %to calculate spike related parameters
        
        
        %Amplitude and time of peaks:
        [peakAmplitude, peakTime] = findpeaks(trace_data(datapulseStart:datapulseEnd-1),'MinPeakDistance',10,'MinPeakProminence',25);
        peaksVector = [peakAmplitude, peakTime]; %vector with the peaks found
        
        datapeakAmplitude{i} = peakAmplitude;
        peakTime = datapulseStart-1 + peakTime; %when adding the peakTime in points to datapulseStart everything shitfs +1 in x, so -1 to shift it back to the right position
        datapeakTime{i} = peakTime;
        
        
        if ~isempty(peakAmplitude)
        derivdata = diff(trace_data,1); %1st derivative to calculate threshold of each AP
        limit = 0.05*max(derivdata); %0.07*max(derivdata); %7% of the AP peak for the threshold to use when 2nd derivative fails
        limit2 = 0.05*max(derivdata); %0.05
        
        troughTime = [];
        troughAmplitude = [];
        thresholdTime = [];
        thresholdTimeSubseq = [];
        thresholdAmplitude = [];
        
        if size(peaksVector,1) > 1
            %Trough:
            %all except the last trough, since it could be after the
            %end of the puls:
            for j = 1:length(peakTime)-1
                [troughAmplitude(j), troughTime(j)] = min(trace_data(peakTime(j):peakTime(j+1)));
                troughTime(j) = peakTime(j)-1 + troughTime(j);
            end
            %now I calculate the last trough
            [troughAmplitudeLast, troughTimeLast] = min(trace_data(peakTime(end):datapulseEnd-1));
            troughTimeLast = peakTime(end) +troughTimeLast;
            %...and add it to the array only if it's before the end of
            %the pulse
            if troughTimeLast < datapulseEnd
                troughAmplitude = [troughAmplitude troughAmplitudeLast];
                troughTime = [troughTime troughTimeLast];
            end
            
            datatroughTime{i} = troughTime;
            datatroughAmplitude{i} = troughAmplitude;
            
            %Onset:
            %first onset:
            %-OPTION 1-
            %using the same method than for 'subsequent APs':
            %[~, thresholdTimeFirst] = min(abs(derivdata(datapulseStart:peakTime(1))-limit));
            %thresholdTimeFirst = datapulseStart-1 + thresholdTimeFirst;
            %-END OPTION 1-
            
            %-OPTION 2-usually works better, because of the capacitive charging of the membrane
            %calculating threshold as the point in which 2nd derivative=0
            %(inflection point):
            A = diff(trace_data(datapulseStart:peakTime(1)),2);
            B = find(A==0); %number of points after datapulse Start where that happens
            
            %If B is empty (no inflection point) then use OPTION 1 code
            if isempty(B) == 0 && length(B) >= 3 % isempty() is not really necessary
                thresholdTimeFirst = datapulseStart-1 + B(end);
            else
                [~, thresholdTimeFirst] = min(abs(derivdata(datapulseStart:peakTime)-limit*1.8)); %1.8 is empirical, otherwise it peaks the threshold too close to peak
                %thresholdTimeFirst = round((peakTime(1)-datapulseStart)/2); %a raw alternative to the 1.8*limit option
                thresholdTimeFirst = datapulseStart-1 + thresholdTimeFirst;
            end
            
            
            %-END OPTION 2-
            
            %subsequent onsets:
            for j = 1:length(peakTime)
                
                if j < length(peakTime) %to avoid going past the last peak
                    [~, thresholdTimeSubseq(j)] = min(abs(derivdata(troughTime(j):peakTime(j+1))-limit2));
                    thresholdTimeSubseq(j) = troughTime(j)-100 + thresholdTimeSubseq(j); %troughTime(j)-10 + thresholdTimeSubseq(j);
                end
            end
            
            thresholdTime = [thresholdTimeFirst thresholdTimeSubseq-10];
            thresholdAmplitude = trace_data(thresholdTime);
            
            datathresholdAmplitude{i} = thresholdAmplitude; %trace_data(thresholdTimeFirst); %
            datathresholdTime{i} = thresholdTime; %thresholdTimeFirst; %thresholdTime;
            
        else %if only one AP
            %Trough:
            [troughAmplitude, troughTime] = min(trace_data(peakTime:datapulseEnd-1));
            troughTime = peakTime-1+troughTime;
            
            datatroughTime{i} = troughTime;
            datatroughAmplitude{i} = troughAmplitude;
            
            %             %Onset:
            %             %-OPTION 1-
            %             %using the same method than for 'subsequent APs':
            %             [~, thresholdTime] = min(abs(derivdata(datapulseStart:peakTime)-limit));
            %             thresholdTime = datapulseStart-1 + thresholdTime;
            %             %-END OPTION 1-
            
            %-OPTION 2-usually works better, because of the capacitive charging of the membrane
            %calculating threshold as the point in which 2nd derivative=0
            %(inflection point):
            A = diff(trace_data(datapulseStart:peakTime(1)),2);
            B = find(A==0); %number of points after datapulse Start where that happens
            
            %If B is empty (no inflection point) then use OPTION 1 code
            if isempty(B) == 0 && length(B) >= 3 % isempty() is not really necessary
                thresholdTime = datapulseStart-1 + B(end);
            else
                %[~, thresholdTime] = min(abs(derivdata(datapulseStart:peakTime)-limit*1.8));
                thresholdTime = round((peakTime(1)-datapulseStart)/2); %a raw alternative to the 1.8*limit option
                thresholdTime = datapulseStart-1 + thresholdTime;
            end
            
            %-END OPTION 2-
            
            thresholdAmplitude = trace_data(thresholdTime);
            
            datathresholdAmplitude{i} = thresholdAmplitude;
            datathresholdTime{i} = thresholdTime;
        end
        
        %-
        %Calculation of avg AP height (peak to trough):
        if size(peaksVector,1) > 1
            AP_height = [];
            ISI = [];
            half_width = []; %these three are defined as empty arrays to ensure that they reset in every cycle, since their size varies!
            %adapt_index = []; %this one doesn't change it's size, and I define it as NaN whenever can't be calculated (when 1 AP only)
            halfTime1 = [];
            halfTime2 = [];
            %AP height:
            if length(troughTime) == length(peakTime) %if it was possible to calculat a trough for every peak (that is, if last trough doesn't fall after pulseEnd)
                for j = 1:length(peakTime)
                    AP_height(j) = abs(peakAmplitude(j)-troughAmplitude(j));
                end
            else
                for j = 1:length(peakTime)-1
                    AP_height(j) = abs(peakAmplitude(j)-troughAmplitude(j));
                end
            end
            
            %ISI:
            for j = 1:length(thresholdTime)
                if j < length(thresholdTime)
                    ISI(j) = thresholdTime(j+1)-thresholdTime(j); %in points
                end
            end
            if length(ISI) > 1
                for k = 1:length(ISI)-1 %this will calculate the differences between ISIs across an AP
                    ISIindex(k) = (ISI(k+1)-ISI(k))/(ISI(k+1)+ISI(k)); %I don't care if points or ms (for which I'd have to do '*dx' above and below, and they'd cancel each other out)
                end
                
                adapt_index = (1/(length(ISI)-1)*sum(ISIindex)); %the closest to 0 the less adaptation in ISIs
            else
                adapt_index = 0; %only applies if more than 2 AP
            end
            %Half-width:
            if length(troughTime) == length(peakTime)
                for j = 1:length(peakTime)
                    t1_data(j) = 0.5*(abs(peakAmplitude(j)-thresholdAmplitude(j)));
                    t1_data(j) = t1_data(j) + thresholdAmplitude(j);
                    [~, halfTime1(j)] = min(abs(trace_data(thresholdTime(j):peakTime(j))-t1_data(j)));
                    halfTime1(j) = thresholdTime(j)-1 + halfTime1(j);
                    
                    t2_data(j) = peakAmplitude(j)-(peakAmplitude(j)-t1_data(j));
                    [~, halfTime2(j)] = min(abs(trace_data(peakTime(j):troughTime(j))-t2_data(j)));
                    halfTime2(j) = peakTime(j)-1 + halfTime2(j);
                    
                    half_width(j) = halfTime2(j)-halfTime1(j);
                end
                
            else
                for j = 1:length(peakTime)-1
                    t1_data(j) = 0.5*(abs(peakAmplitude(j)-thresholdAmplitude(j)));
                    t1_data(j) = t1_data(j) + thresholdAmplitude(j);
                    [~, halfTime1(j)] = min(abs(trace_data(thresholdTime(j):peakTime(j))-t1_data(j)));
                    halfTime1(j) = thresholdTime(j)-1 + halfTime1(j);
                    
                    t2_data(j) = peakAmplitude(j)-(peakAmplitude(j)-t1_data(j));
                    [~, halfTime2(j)] = min(abs(trace_data(peakTime(j):troughTime(j))-t2_data(j)));
                    halfTime2(j) = peakTime(j)-1 + halfTime2(j);
                    
                    half_width(j) = halfTime2(j)-halfTime1(j);
                end
                
            end
            %halftime1(i) = halfTime1;
            %halftime2(i) = halfTime2;
        else
            %AP height and ISI:
            if length(troughTime) == length(peakTime)
                AP_height = abs(peakAmplitude-troughAmplitude);
            else
                AP_height = NaN;
            end
            ISI = NaN;
            adapt_index = NaN;
            
            %Half-width:
            t1_data = 0.5*(abs(peakAmplitude-thresholdAmplitude));
            t1_data = t1_data + thresholdAmplitude; %target mV to look for t1
            
            [~, halfTime1] = min(abs(trace_data(thresholdTime:peakTime)-t1_data));
            halfTime1 = thresholdTime-1 + halfTime1;
            
            t2_data = peakAmplitude-(peakAmplitude-t1_data); %same as t1_data in fact, but nice to have it in case I want to plot it
            [~, halfTime2] = min(abs(trace_data(peakTime:troughTime)-t2_data));
            halfTime2 = peakTime-1 + halfTime2;
            
            half_width = halfTime2-halfTime1; %in points
            
        end
        dataAPheight{i} = AP_height;
        dataISI{i} = ISI;
        datahalfwidth{i} = half_width;
        
        
        datahalftime1{i} = halfTime1;
        datahalftime2{i} = halfTime2;
        
        
        spike_result(i,3) = NaN;
        spike_result(i,4) = size(peaksVector,1); %stores number of detected APs in each trace
        spike_result(i,5) = size(peaksVector,1)/(length(datapulseStart:datapulseEnd-1)*dx/1000); %spike freq in spikes/sec for f-i curve; denominator = 1 sec
        spike_result(i,6) = datathresholdAmplitude{i}(1); %onset amplitude (mV) for 1st spike
        spike_result(i,7) = (datathresholdTime{i}(1)-datapulseStart)*dx; %onset latency (ms) for 1st spike
        spike_result(i,8) = datapeakAmplitude{i}(1); %absolute peak amplitude (mV) for 1st spike
        spike_result(i,9) = (datapeakTime{i}(1)-datapulseStart)*dx; %peak latency (ms) for 1st spike
        spike_result(i,10) = dataAPheight{i}(1); %AP height (mV) for 1st AP in the step
        spike_result(i,11) = dataAPheight{i}(end); %AP height (mV) for the last AP, if there's more than one AP in the step
        spike_result(i,12) = mean(dataAPheight{i}); %avg AP height (mV)
        spike_result(i,13) = dataISI{i}(1)*dx; %1st ISI (between 1st two APs) in ms
        spike_result(i,14) = mean(dataISI{i}*dx); %avg ISI in ms
        spike_result(i,15) = adapt_index; %it tells us something about how ISI changes from AP to AP (is the train accelerating or not?)
        spike_result(i,16) = mean(datahalfwidth{i}*dx); %avg HW in ms
        
        end 
    end 
    
end

spike_result_str = {'I_step', 'Vm', 'SAG_mV', 'Spike_count', 'Spikes/sec','ThresholdAmp_1st_AP_mV','ThresholdLatency_1stAP_ms','PeakAmp_1st_AP_mV','PeakLatency_1stAP_ms', 'AP_height_1stAP', 'AP_height_lastAP', 'avgAP_height_mV','First_ISI_ms','avgISI_ms','Adapt_index','avgHW_ms'};
%if there were real data, convert the cell array into a mat to be stored in spikes table
spikes_table = array2table(spike_result,'VariableNames',spike_result_str);

    % %%
% %Plot all the traces together:
% for i = 1:size(data,2)
%
%     subplot1 = subplot(5,10,i);
%     trace_data = smooth(data(:,i));
%
%     if i >= obs_threshold
%         plot(time,trace_data,'k','LineWidth',0.5,'Parent',subplot1);
%         hold on; plot(datapeakTime{i}*dx, datapeakAmplitude{i}, 'bo','Parent',subplot1)
%         hold on; plot(datatroughTime{i}*dx, datatroughAmplitude{i}, 'bo','Parent',subplot1)
%         hold on; plot(datathresholdTime{i}*dx, datathresholdAmplitude{i}, 'ro','Parent',subplot1)
%         hold on; plot(datahalftime1{i}*dx, trace_data(datahalftime1{i}),'go','Parent',subplot1)
%         hold on; plot(datahalftime2{i}*dx, trace_data(datahalftime2{i}),'go','Parent',subplot1)
% %         hold on;
% %         plot((dataRCpulseStart:dataRCpulseStart+80/dx)*dx,realfitfallingPhase_all{i}+baseline(i),'g','LineWidth',2)
%     elseif i < obs_threshold
%         plot(time,trace_data,'k','LineWidth',0.5,'Parent',subplot1);
%         hold on;
%         if i < indexZero
%             plot(dataminTime(i)*dx, dataminAmplitude(i), 'ro', 'Parent',subplot1);
%             plot(datasteadyTime(i)*dx, datasteadyAmplitude(i), 'bo', 'Parent',subplot1);
%         end
% %         hold on;
% %         plot((dataRCpulseStart:dataRCpulseStart+80/dx)*dx,realfitfallingPhase_all{i}+baseline(i),'g','LineWidth',2)
%
%     end
% xlabel('Time (msec)');
% ylabel('Membrane potential (mV)'); hold on;
% maxvalue = max(trace_data(4000:end));
% minvalue = min(trace_data(4000:end));
% ylim([minvalue-50 maxvalue+50]);
% end
%%
%Plot traces individually:

for i= 1:size(dataTraces,2)%obs_threshold:obs_threshold+1
    trace_data = smooth(dataTraces(:,i));
    if i < indexZero
        figure;
        plot(time,trace_data,'k','LineWidth',0.5);
        hold on;
        plot(dataminTime(i)*dx, dataminAmplitude(i), 'ro');
        plot(datasteadyTime(i)*dx, datasteadyAmplitude(i), 'bo');
    elseif i >= rheoIndex
        figure;
        plot(time,trace_data,'k','LineWidth',0.5);
        hold on;
        plot(datapeakTime{i}*dx, datapeakAmplitude{i}, 'bo')
        plot(datatroughTime{i}*dx, datatroughAmplitude{i}, 'bo')
        plot(datathresholdTime{i}*dx, datathresholdAmplitude{i}, 'ro')
        plot(datahalftime1{i}*dx, trace_data(datahalftime1{i}),'go');
        plot(datahalftime2{i}*dx, trace_data(datahalftime2{i}),'go');
        %      hold on;
        %          plot((dataRCpulseStart:dataRCpulseStart+80/dx)*dx,realfitfallingPhase_all{i}+baseline(i),'g','LineWidth',2)
    end
    xlabel('Time (msec)');
    ylabel('Membrane potential (mV)'); hold on;
    maxvalue = max(trace_data(4000:end));
    minvalue = min(trace_data(4000:end));
    ylim([minvalue-50 maxvalue+50]);
    
    if plotPause
        pause;
    end 
    close all
end
%% store results
save(fullfile(pathname,[extractBefore(fileName, '.txt') '_SpikeTable.mat']),'spikes_table');

end


%% lucas Store results:
%resultvar = [pathname fileName(1:end-4) '_result.mat'];
%save(resultvar,'IV_results');

% argh = [fileName(1:end-4) '_result.mat'];
% resultvar = [pathname argh];
% save(resultvar, results)

% %...and print these:
% avg_Vm = mean(IV_results.Vm, 'omitnan');
% disp(['mean Vm for all traces is: ' num2str(avg_Vm) 'mV']);
% disp(['rheobase is: ' num2str(threshold_current) 'pA']);
% disp(['Input resistance is: ' num2str(Rm) 'MOhms']);
% disp(['Series resistance is: ' num2str(Rs) 'MOhms']);


