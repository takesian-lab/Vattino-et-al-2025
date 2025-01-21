function [projData] = writeProtInProjData(meta,cellID,protNum,resTableName,resTable,varargin)
% input parser for opting to using preloaded projData or newly loaded projdata  
p = inputParser; 
%USAGE: addOptional(p,'parametername',defaultvalue);
% set default values for the possible input variables
addOptional(p, 'inputProjData', ''); % this default is loaded projData based on meta data

parse(p, varargin{:});
ops = p.Results; 
% ------- end parse varargin

if ~isempty(ops.inputProjData)
    projData = ops.inputProjData;
else 
    load(meta.save_file);% load projData from the metafile (file location is stored in meta.save_file)
end 


%% load dataa and find cooresponding cell 
EphysDate = meta.EphysDate;

indexOG = find((ismember({projData.EphysDate}, EphysDate) & ismember({projData.CellID}, cellID))==1);
if length(indexOG) == 1
    index = indexOG;
else
    error('STOP: more than one corresponding cell found in data structure')
end

%% write into projData struct
% if there are no such field, write the current data directly into a new field in this struct
if ~isfield(projData,[resTableName])
    %projData(index).minStimThreshold = stimThreshold;
    projData(index).(resTableName) = resTable;
    
    % if there are fields for Results, then test if there is already a table in this field
elseif ~istable(projData(index).(resTableName)) %
    %projData(index).minStimThreshold = stimThreshold;
    projData(index).(resTableName) = resTable;
    
    % if there are no duplicate prot num, then add onto it
elseif  ~any(ismember(projData(index).(resTableName).protNum,protNum)) % and if there's NO duplicate protocal number
    projData(index).(resTableName) = [projData(index).(resTableName);resTable]; %add the current baselineTable to the existing table
    
else %if there are duplicates protocols, then replace it with the current one
    duplicate=ismember(projData(index).(resTableName).protNum,protNum);
    disp('Duplicate found in projData, replacing it with the current analysis...');
    projData(index).(resTableName)(duplicate,:) = resTable;
end
% save(save_file,'projData','-v7.3');

end 