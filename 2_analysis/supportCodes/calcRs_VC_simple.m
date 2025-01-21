function [Rin, Rs, Cm, error] = calcRs_VC(data, acq_freq, stimTable)
% function for calculating Rs for voltage clamp experiments
% Input: 
%   1. data: single data trace
%   2. acq+freq: acquisition frequence (dx), 
%   3. stimTable: RC stim settings for a single input -- could be already loaded in the prot folder with data or read from the raw data trace's user data (need to be loaded before)  a
% Output:
%   1. Rin
%   2. Rs
%   3. Cm
%   4. error
%% define stim variables from the stim table for this acq
pulseAmp = stimTable.amplitude;
pulseStart = stimTable.delay;
pulseStop = stimTable.delay + stimTable.pulseWidth; 
baselineStart = 1;
baselineEnd = stimTable.delay-1;

%% define stim variables from the stim table for this acq
pulseAmp = stimTable.amplitude;
pulseStart = stimTable.delay;
pulseStop = stimTable.delay + stimTable.pulseWidth; 
baselineStart = 1;
baselineEnd = stimTable.delay-1;
Rin=0;
Rs=0;
Cm=0;
error=0; %will turn into 1 if rin, rs and cm can't be calculated

pulseLength = pulseStop-pulseStart;

% Time variables in sampling points:
databaselineStart = round(baselineStart/acq_freq); 
databaselineEnd = round(baselineEnd/acq_freq);
datapulseStart = round(pulseStart/acq_freq); 
datapulseStart = datapulseStart+1;
datapulseStop = round(pulseStop/acq_freq);
datapulseLength = round(pulseLength/acq_freq);


if pulseAmp>0
    [peak, peakPos] = max(data(datapulseStart:datapulseStop));
else
    [peak, peakPos] = min(data(datapulseStart:datapulseStop));
end
peakPos = peakPos + datapulseStart-1;

%delta=round(0.1*pulseLength/dx);
%endline corresponds to the I in the ohmic phase of the V drop:
endline = mean(data(datapulseStart + round(datapulseLength*0.7):datapulseStop-1)); 
%baseline is pre-pulse current:
baseline = mean(data(databaselineStart+1:databaselineEnd-1));

%Looking for data points above endline + 1/5 of endline-to-peak difference:
if pulseAmp>0
    above=find(data(datapulseStart:datapulseStop)>((peak-endline)/5)+endline);
    if isempty(above)
        error=1;
        return;
    end
    last=above(end)+datapulseStart;
else
    above=find(data(datapulseStart:datapulseStop)<((peak-endline)/5)+endline);
    if isempty(above)
        error=1;
        return;
    end
    last=above(end)+datapulseStart;
end
delta=round((last-peakPos)/2); 
%delta in points between last point that was above endline + 1/5 of
%endline-to-peak difference and the actual peak position
clear above last

peak1=peak-endline;
peak2=data(delta+peakPos)-endline;
peak3=data(2*delta+peakPos)-endline;
peakPos=peakPos-datapulseStart-1;

if (peak1-peak2)*(peak2-peak3)<0
    disp('calcRs: 3 points are non-monotonic.  No fit possible')
    error=1;  
    Rin = NaN;
    Rs = NaN;
    Cm = NaN;
    return
end

tau=delta*(1/log(peak1/peak2)+1/log(peak2/peak3)+2/log(peak1/peak3))/3;
amp=(peak1/exp(-peakPos/tau)+peak2/exp(-(peakPos+delta)/tau)+peak3/exp(-(peakPos+2*delta)/tau))/3;

Rs=1000*pulseAmp/amp;
Rin=1000*pulseAmp/(endline-baseline)-Rs;
Cm=1000*tau*acq_freq/Rs;


if any(~isreal([tau amp Rs Rin Cm]))
    disp('calcRs: complex number returned.  No fit possible')
    Rin=NaN;
    Rs=NaN;
    Cm=NaN;
    error=1;
    return
end



clear tau amp peak peak1 peak2 peak3 peakloc endline baseline
end