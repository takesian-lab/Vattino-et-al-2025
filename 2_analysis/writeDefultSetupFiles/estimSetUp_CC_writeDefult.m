%% 4. Analysis -  establish opto settings
%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
clear all
dx = 0.1;

%Parameters for the stimulation protocol (in ms):
searchWindow = 50; %search peaks in the x msec  following stim start
searchAfterStim = 2; %when to start searching, after the stimStart in msec (3 sec for estim) 
timePreStim = 50; %time before stim that we will use to average baseline
RsThreshold = 50; %series resistant threshold: above which that trace will be considered unhealth (if each trace recorded RC check)

%% Peak detection
minPeakDist = 1.5./dx; % distance between two synaptic events in msec (converted to datapoints)
minPeakProm = 0.1; % in mV (how big of an amplitude difference in mV for two "peaks" - 0.1 for disynaptic EPSPs in christine's MGv stim protcols

% how many peaks are you looking to detect (synaptic events) for single
% stimulation 
numPeaks = 2;

%Criteria to detect events:
% max msec after stim(converted to data points) for first and second peaks
% (second peak not required if numPeaks = 1;
onsetTime_ref_peak1 = 7./dx; %5 for Christine's MGv monosynaptic EPSPs/APs
onsetTime_ref_peak2 = 30./dx; %10 for Christine's MGv disynaptic EPSPs

limitAPthresh = 0.05; % porpotion of the max AP height use as the AP threshold if derivative fails
onsetCheckAmp = 0.1; %%look for values in the data that are less than 5pA (or 0.01mV) away from ideal 20%



%% old ones 
peakAmplitudeRef = 1; %range 5-8 pA
% for estimating onset 
lowEst = 0.3; %to estimate the onset of mono peak, draw linear line from 20% amp
highEst = 0.8; %to estimate the onset of mono peak, draw linear line to 70% amp

% save the settingFile
%settingFile = fullfile('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings','monoSetup_default.mat');
save('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\eStimSetup_CC_default.mat');