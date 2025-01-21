function [allStimTable,allStimTraces] =  loadStim(protPath)
%loadStim: 
% input:  pathname (str): could be pre-defined in meta.pathname OR
% generated with a function - [pathname] = findProtPath(cellID,protNum,EphysDate,base_dataPath)
% output: 1. allStimTable (a structure with fields ao0,ao1,ao2,ao3: each field is a table that stores all the variables of the stim params)
%         2. allStimTraces (a structure with fields ao0,ao1,ao2,ao3: each
%         field is a double array)

%% Load data file to analyze
file = dir([protPath '/*.txt' ]); % to find the text file with data

% if it has already been analyzed (already have a second txt file or have other stim.txt type files), use the first one, which is the original raw data set
sortedFileName = sort({file.name});
dataFileName = sortedFileName{1}; %the first .txt file should have nothing added to Acqx_Acqy

%% load stim struct (or create if not already exist)
pulseStructFile = dir([protPath '/*_PulseParamsStruct.mat' ]);
pulseTraceFile =  dir([protPath '/*_PulseTraces.mat' ]);
%
if isempty(pulseStructFile) %if no pulse struct already exist in this folder, create one   
    epoacqrange_temp = regexp(dataFileName,'\d*','Match');
    epoacqrange = [str2num(epoacqrange_temp{1}):str2num(epoacqrange_temp{2})];
    
    %use deconstructProtPath to find the useful info from the pathname

    [~,~,cellID,protName,~]=deconstructProtPath([protPath '\']); 
    ephysDayPath = extractBefore(protPath,cellID);
    dataPath = extractBefore(protPath,protName);
    
    %use extractStimPattern function to read and store allStimTable and
    %allStimTraces
    [allStimTable,allStimTraces] = extractStimPattern(ephysDayPath,dataPath,protName,epoacqrange);
else    
    %load allStimTable
    load([pulseStructFile.folder '\' pulseStructFile.name]); %this will load the file, which will appear as 'allStimTable' (struct)
    %load pulseTraceFile
    load([pulseStructFile.folder '\' pulseTraceFile.name]);
    
end    
%% compare the fieldnames of the allStimTable
if any(contains(fieldnames(allStimTable),'ao')) 
    allStimTableTemp = allStimTable;
    clear allStimTable
    % load meta data from raw data folder
        metafolder = fileparts(fileparts(pulseStructFile.folder));
        fileListing = dir(metafolder);
        metaFile =  fileListing(find(contains({fileListing.name},'metaData')));
        load(fullfile(metaFile.folder,metaFile.name));
        
        % check if there's already aoMatch filed in meta. If not create
        % default and save the edited meta struct
        if ~isfield(meta,'aoMatch')
            meta.aoMatch.ao0 = 'intrinsic';
            meta.aoMatch.ao1 = 'led';
            meta.aoMatch.ao2 = 'estim';
            meta.aoMatch.ao3 = '';
            save(meta.metaFilePath, 'meta');
        end
        
        % now rewrite the name
        aoFields = fieldnames(meta.aoMatch);
        for f = 1:length(aoFields)
            if ~isempty(meta.aoMatch.(aoFields{f}))
                allStimTable.(meta.aoMatch.(aoFields{f})) = allStimTableTemp.(aoFields{f});
            end
        end
        allStimTable.RC = allStimTableTemp.RC;
end

%% compare the fieldnames of the allStimTrace
if any(contains(fieldnames(allStimTraces),'ao')) 
    allStimTracesTemp = allStimTraces;
    clear allStimTraces
        % now rewrite the name
        aoFields = fieldnames(meta.aoMatch);
        for f = 1:length(aoFields)
            if ~isempty(meta.aoMatch.(aoFields{f}))
                allStimTraces.(meta.aoMatch.(aoFields{f})) = allStimTracesTemp.(aoFields{f});
            end
        end
end
end 
