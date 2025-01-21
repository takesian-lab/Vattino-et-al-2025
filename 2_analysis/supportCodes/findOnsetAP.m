function [onsetAmp,onsetPoint] = findOnsetAP(data,datasearchStart,dataThisPeakTime,limitAPthresh)
% find onsetAmp and onsetTime (number of data points from the start of trace) if this is the first AP in the trace. 

A = diff(data(datasearchStart:dataThisPeakTime),2); % find the second derivvative of data
B = find(A==0); %find inflection point: number of data points after datapulse Start where that happens

derivdata = diff(data,1); %1st derivative to calculate threshold of each AP

%If B is empty (no inflection point) then
if isempty(B) == 0 && length(B) >= 3 
    onsetTime_temp = B(end)-1; %number of data points since search start
else
    [onsetAmp, onsetTime_temp] = min(abs(derivdata(datasearchStart:dataThisPeakTime)-limitAPthresh*1.8)); %1.8 is empirical, otherwise it peaks the threshold too close to peak
    onsetTime_temp = onsetTime_temp-1; %number of data points since search start
end

% convert to data point since the beginning of the trace
onsetPoint = onsetTime_temp+datasearchStart;
