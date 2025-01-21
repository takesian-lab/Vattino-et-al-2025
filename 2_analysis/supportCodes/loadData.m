function [traces,fileName] =  loadData(pathname)
%loadData: 
% input:  pathname (str): could be pre-defined in meta.pathname OR
% generated with a function - [pathname] = findProtPath(cellID,protNum,EphysDate,base_dataPath)
% output: traces (double array)

%% Load data file to analyze
A = dir(pathname); 
file = A(contains({A.name},'.txt'));
% file = dir([pathname '/*.txt' ]); % to find the text file with data

% if it has already been analyzed (already have a second txt file or have other stim.txt type files), use the first one, which is the original raw data set
sortedFileName = sort({file.name});
fileName = sortedFileName{1}; %the first .txt file should have nothing added to Acqx_Acqy 

traces = load(fullfile(pathname,fileName)); % load traces
end 