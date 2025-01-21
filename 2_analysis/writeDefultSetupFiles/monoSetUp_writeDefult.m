%% 4. Analysis -  establish monosynaptic settings
%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
clear all
dx = 0.1;

%Parameters for the stimulation protocol (in ms):
searchWindow = 15; %seach peaks in the x msec following stim start... for both mono and disynaptic peaks
searchAfterStim = 5; %when to start searching, after the start of the stim artifact (artifact should be 0.5 ms long)
timePreStim = 50; %time before stim that we will use to average baseline (this means we wil look from 

%Criteria to detect events - could be modified 
% stimulation protocols:
peakAmplitudeRef = 7; %range is 5-8 pAs
doublePeakDist = 40; %max distance between two peaks for mono and disynaptic peaks

% for estimating onset mono 
lowEst = 0.25; %to estimate the onset of mono peak, draw linear line from 25% amp
highEst = 0.75; %to estimate the onset of mono peak, draw linear line to 25% amp

RsThreshold = 100; % the series resistance threshold to have the epoch thrown out

% save the settingFile
%settingFile = fullfile('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings','monoSetup_default.mat');
save('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\monoSetup_default.mat');
