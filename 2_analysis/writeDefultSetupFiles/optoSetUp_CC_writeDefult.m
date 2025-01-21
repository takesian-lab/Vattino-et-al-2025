%% 4. Analysis -  establish opto settings
%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
clear all
dx = 0.1;

%Parameters for the stimulation protocol (in ms):
searchWindow = 30; %search peaks in the x msec following stim start
searchAfterStim = 1; %when to start searching, after the stimStart
timePreStim = 50; %time before stim that we will use to average baseline (this means we wil look from 
RsThreshold = 100; % the series resistance threshold to have the epoch thrown out

%Criteria to detect events - could be modified 
% stimulation protocols:
peakAmplitudeRef = 8;
onsetTime_ref_EPSP = 30;
onsetTime_ref_IPSP = 40;
fallTimefromPeak_ref = 3; %in ms
totalCharge_ref_IPSP = 0.01; %ref charge for IPSCs
totalCharge_ref_EPSP = -0.09; %ref charge for EPSCs - previously -0.09
fallTimeMax = 2000; %max possible return to baseline - previously 1200
% for estimating onset 
lowEst = 0.3; %to estimate the onset of mono peak, draw linear line from 20% amp
highEst = 0.8; %to estimate the onset of mono peak, draw linear line to 70% amp

% save the settingFile
%settingFile = fullfile('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings','monoSetup_default.mat');
save('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\optoSetup_CC_default.mat');