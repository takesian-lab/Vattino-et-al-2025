%% setup params for train 

% TODO: need to be changed for train data

%Sampling rate in ms: 0.1 ms = 100 us (10kHz - 10 points per ms)
clear all
dx = 0.1;

%Parameters for the stimulation protocol (in ms):
searchWindow = 30; %search peaks in the x msec following stim start
searchAfterStim = 1; %when to start searching, after the stimStart
timePreStim = 50; %time before stim that we will use to average baseline (this means we wil look from 
RsThreshold = 200;
%Criteria to detect events - could be modified 
% stimulation protocols:
startSearch = 8*1/dx; %2.5 ms or 25 points after repsonse; 8 ms or 80 points after response
stopSearch = 30*1/dx;

% save the settingFile
%settingFile = fullfile('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings','monoSetup_default.mat');
save('\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\trainSetup_default.mat');