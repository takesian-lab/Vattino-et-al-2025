function [onsetAmp, onsetPoint] = findOnsetTime(data,dataThisPeakTime,stimStartPoint,Setup)
% find onsetAmp (real voltage) and onsetPoint (# data points since the very first data point) with the walkback method for any subthreshold (not AP)
% calculations

% write default values from setup

dx = Setup.dx;
datasearchStart = stimStartPoint + Setup.searchAfterStim./dx;
lowEst = Setup.lowEst;
highEst = Setup.highEst;
onsetCheckAmp = Setup.onsetCheckAmp;
stimStart = stimStartPoint;

% set default onsetAmp and onsetTime 
onsetAmp = data(datasearchStart);
onsetPoint = datasearchStart;
%% find 20% estimate
ideal_twenty_percent = lowEst*(data(dataThisPeakTime)-data(datasearchStart)); %real 20% value
twenty_percent_amp = ideal_twenty_percent+data(datasearchStart); %the amplitude of that ideal 20%
k = dataThisPeakTime;
k2=k-length(datasearchStart:dataThisPeakTime)+1; %k-length(data_amp);
foundk = 0;

while k > k2
    if abs(twenty_percent_amp-data(k)) > onsetCheckAmp %look for values in the data that are less than 5 pA away from ideal 20%; 5 pA can be modified
        k = k-1;
        foundk = 0;
    else
        foundk = 1;
        %k=k;
        break
    end
    
end
twenty_percent_amp = data(k);
twenty_percent_time = k*dx;


%to find 70% value in data: walk forward from the 20% time until the peak time
ideal_seventy_percent = highEst*(data(dataThisPeakTime)-data(datasearchStart)); %real 20% value
seventy_percent_amp = ideal_seventy_percent+data(datasearchStart); %the amplitude of that ideal 20%
j = k;%twenty_percent_time(i)/dx;
j2 = j+length(twenty_percent_time/dx:dataThisPeakTime)-1;
foundj = 0;
while j < j2 %j <= j2
    if abs(seventy_percent_amp-data) > 8 %look for values in the data that are less than 5 pA away from ideal 20%; 5 pA can be modified
        j = j+1;
        foundj = 0;
    else
        foundj = 1;
        %j=j;
        break
    end
end

seventy_percent_amp = data(j);
seventy_percent_time = j*dx;

x1 = k:round(dataThisPeakTime); %or finalFirstPeakTime(i)/dx
%     x1 = k:j;
if length(x1) >= 2 %requirement of 2 points or more to fit rise; if id doesn't find points everything crashes (this happens when responses are abolished)
    
    FR = fit(x1', data(x1), 'poly1'); %same as using polyfit(x1',data(x1),1) --> linear fit
    CR = coeffvalues(FR);
    %lin_FR_vec = k:realpeakTime(i);
    lin_FR_vec = k:j;
    lin_FR_dist = 100-length(lin_FR_vec); %I want 100 pts lines, so this is the number of points to add
    %lin_FR = CR(1)*(k-lin_FR_dist:realpeakTime(i))+CR(2); %this has 100 pts
    lin_FR = CR(1)*(k-lin_FR_dist:j)+CR(2); %this has 100 pts
    lin_FR_plot = lin_FR;
    %lin_FR_time = realpeakTime(i)-length(lin_FR)+1:realpeakTime(i); %and this is the time axis
    lin_FR_time = j-length(lin_FR)+1:j; %and this is the time axis
    lin_FR_time_plot = lin_FR_time;
    x2 = stimStart-2000:stimStart-10;
    FB = fit(x2', data(x2),'poly1');
    CB = coeffvalues(FB);
    lin_FB = CB(1)*(datasearchStart-5:datasearchStart+length(lin_FR)-6)+CB(2);
    lin_FB_plot = lin_FB;
    lin_FB_time = datasearchStart-5:datasearchStart+length(lin_FB)-6;
    lin_FB_time_plot = lin_FB_time;
    
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
        onsetAmp = data(k);
        plot_onset_point = k*dx; %added (LV_01132021)
        onsetPoint=plot_onset_point/dx; %added (LV_01132021)
       
        
    else
        
        for k = 1:length(t_vector)
            closest_to_intersect(k) = min(abs(fit_onset_amp-data(t_vector(k))));
        end
        
        %among values in closest_to_intersect I want theclosest to 0 with highest index:
        [min_val,min_indices] = mink(closest_to_intersect,3);
        onset_point = t_vector(max(min_indices)); %point time in which the onset happens
        
        onsetAmp = data(onset_point);
        onsetPoint= onset_point; %real onset latency
        
    end
    
end
% figure;
% plot(data)
% hold on
% scatter(twenty_percent_time./dx,twenty_percent_amp,'k'); hold on
% scatter(seventy_percent_time./dx,seventy_percent_amp); hold on
% scatter(onset_point,onsetAmp);
end