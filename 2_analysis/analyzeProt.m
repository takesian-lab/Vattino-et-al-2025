function [projData] = analyzeProt(protPath,varargin)
% VC analyis:
% 1: minimal e stimulation analysis (detectSinglePeak.m) - no charge calc
% 2: single stim analysis (either estim or led)
% 3: train stim (either estim or opto)
% 4: train ramp (Christine's thalamic cells)


% CC analysis
% 11: full FI curve (full IV curve or single current injection)
% 12: hyperpolarizing currents over long period of time
% 13: single stim CC


%% load metadata and project data
[~,EphysDate,cellID,protName,protNum]=deconstructProtPath(protPath);

EphysDayPath = extractBefore(protPath,cellID);
cd(EphysDayPath);
metaFileStruct = dir('*metaData.mat');

if isempty(metaFileStruct)
    error('no meta data foulnd in Ephys Day folder. Check to continue');
else
    load(metaFileStruct.name); %load meta
end


%% ------- parse varargin
p = inputParser;
%USAGE: addOptional(p,'parametername',defaultvalue);

% set default values for the possible input variables
addOptional(p, 'plotPause', 0);
addOptional(p, 'reAnalyze', 0);
addOptional(p, 'inputProjData', ''); % this default is loaded projData based on meta data

parse(p, varargin{:});
ops = p.Results;
% ------- end parse varargin

plotPause = ops.plotPause;
reAnalyze = ops.reAnalyze;

if ~isempty(ops.inputProjData)
    projData = ops.inputProjData;
else
    A = load(meta.save_file); %load projData
    projData = A.projData;
end
%% find corresponding protocal from protocal table

% if protList excel file doesn't excist yet, create one
if ~isfield(meta,'protListFile')
    [~,~]=createProtList(meta.pathname,meta.MouseID);
end

protTable = readtable(meta.protListFile);
ind = find(protTable.ProtNum==protNum & contains(protTable.CellID,cellID));

%%

switch(protTable.analysisIndex(ind))
    
    case 1 % 1: VC - Min eStim evoked
        %% 1 = VC - Min eStim evoked
        disp(['Evaluating Min eStim evoked (#1) for ' EphysDate ' ' cellID ' : ' protName ' ...']);
        resTableName = 'minResults';
        if (findProtInProjData(projData,protPath,'minResults') && reAnalyze ) ... % if we can find the same prot alareay analyzed and that we want to reanalyze
                ||  ~findProtInProjData(projData,protPath,'minResults') % or we didn't find it to be alareayd analyzed
            if reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Reanalyzing: ' protName '...']);
            end
            % TODO: add polarity for cellIDs bla and prot number bla --> call SetupFile in excel sheet
            % TODO: TO TEST THIS OUT
            %3.1 load meta file and edit:
            %default file should be same for all users for similar analysis
            monoSetup = load(protTable.SetupFile{ind});
            disp(monoSetup);%display monoSetup to check for any alterations
            disp('check for any changes in command window. To edit: monoSetup.xxx = nnn; When finished, write "dbcont"');
            %
            keyboard
            stimThreshold = input(['minimal stim threshold for "' protName '"is:']);
            polarity = -1;
            
            [resultTableDraft,traces,stimTraces,plottingStruct] = detectMinPeak(cellID,protNum,meta,polarity,monoSetup);
            %this range (in datapoints) is preset as the xlim in the app.
            run('traceByTrace_viewUpdate.mlapp')
            pause;
            [resTableName,resTable]= minStim_updateSave(cellID,protNum,stimThreshold,resultTableDraft,traces,stimTraces,meta,monoSetup);
            
            % To save struct after analysis into projData
            [projData] = writeProtInProjData(meta, cellID,protNum,resTableName,resTable,'inputProjData',projData);
            clear resTable resultTableDraft
        else
            disp(['Already analyzed, not reanalyzing: ' protName '...'])
        end
        
    case 2 % VC - Single stim evoked (Opto or estim)
        %% 2. VC - Single stim evoked (Opto or estim)
        disp(['Evaluating VC single stim evoked (#2) for ' EphysDate ' ' cellID ' : ' protName ' ...']);
        
        stimType = protTable.StimType{ind}; %define stim type -- estim, led488, picosprtizer -- should match the ao outputs
        resTableName = 'singleStimResults';
        
        % --- find if the data already exists in project data ---
        if (findProtInProjData(projData,protPath,'Results') && reAnalyze ) ... % if we can find the same prot alareay analyzed and that we want to reanalyze
                ||  ~findProtInProjData(projData,protPath,'Results') % or we didn't find it to be alareayd analyzed
            if reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Reanalyzing: ' protName '...']);
            elseif ~reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Already analyzed, not reanalyzing: ' protName '...'])
            end
            
            % --- load setup ---
            Setup = load(protTable.SetupFile{ind});
            
            if plotPause % display and allow changes to setup params
                disp(Setup);%display optoSetup to check for any alterations
                disp('check for any changes in command window. To edit: optoSetup.xxx = nnn; When finished, write "dbcont"');
                keyboard
            end
            
            % ------------ write polarity ------------
            %check for Recording voltage in protocal table
            [polarity] = writePolarityFromProtTable(protTable,ind);
            
            % -----------write if TTX_4AP (1 or 0) and change the search window in setup value---
            if contains(stimType,{'led','opto'},'Ignorecase',1) % only check for TTX and 4AP for opto experiments
                if any(strcmp('DrugUsed', protTable.Properties.VariableNames))
                    if  ~ismissing(protTable.DrugUsed(ind)) % if protocal table has the column "DrugUsed"
                        TTX_4AP = contains (protTable.DrugUsed(ind),'TTX') || contains (protTable.DrugUsed(ind),'4AP');
                    else  % if not in protocal table or data entry not numeric, ask for manual input
                        TTX_4AP = input(['added TTX and 4AP for "' protName '? Yes(1) or No(0)...']);
                    end
                end
            end
            % increase the searchwindow if has TTX and 4AP
            if TTX_4AP == 1; Setup.searchWindow = 45; end
            
            % ----------- analyze each trace and visualize each trace with fit --=
            [resultTableDraft,traces,stimTraces] = detectStim_VC(cellID,protNum,meta,polarity,Setup,stimType);
            
            % if user wants to stop and look at all the traces, pause code for user to inspect
            if plotPause; pause;end;   close all; 
            
            % ----------- write custom vars and update all the calculated numbers
            if contains(stimType,'led','Ignorecase',1)
                var1 = protTable.DrugUsed{ind}; var1Name = 'DrugUsed';
                var2 = protTable.StimLevel(ind); var2Name = 'StimLevel';
            elseif contains(stimType,'estim','Ignorecase',1)
                var1 = polarity; var1Name = 'polarity';
                var2 = protTable.StimLevel{ind}; var2Name = 'stimLevel';
            end
            [resTableName,resTable]= stim_updateSave(cellID,protNum,var1,var1Name,var2,var2Name,resTableName,resultTableDraft,traces,stimTraces,meta,Setup);
            
            % ------------To save struct after analysis into projData
            [projData] = writeProtInProjData(meta, cellID,protNum,resTableName,resTable,'inputProjData',projData);
            
        end
        
    case 3 % train stims: either estim (default) or opto stim
        %% 3: train stims: either estim (default) or opto stim
        disp(['Evaluating train (#5) for ' EphysDate ' ' cellID ' : ' protName ' ...']);
        resTableName = 'trainResults';
        if (findProtInProjData(projData,protPath,resTableName) && reAnalyze ) ... % if we can find the same prot alareay analyzed and that we want to reanalyze
                ||  ~findProtInProjData(projData,protPath,resTableName) % or we didn't find it to be alareayd analyzed
            if reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Reanalyzing: ' protName '...']);
            elseif ~reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Already analyzed, not reanalyzing: ' protName '...'])
            end
            trainSetup = load(protTable.SetupFile{ind});
            % --- write polarity ---
            %check for Recording voltage & polarity in protocal table            
            [polarity] = writePolarityFromProtTable(protTable,ind);
          
            % ---- write drug used
            drug = protTable.DrugUsed{ind};
            
            % --- analyze each trace and visualize each trace with fit ---
            stimtypeInd = find(contains(protTable.Properties.VariableNames,'StimType'));
            
            if ~isempty(stimtypeInd)
                stimType = protTable.StimType{ind};
            else
                stimType = 'estim';
            end

            [resTableName,resTable] = detectTrain(cellID,protNum,meta,polarity,trainSetup,stimType,drug);

            % To save struct after analysis into projData
            [projData] = writeProtInProjData(meta, cellID,protNum,resTableName,resTable,'inputProjData',projData);            
        end
        
        if plotPause
            pause
        end
        close all

    case 4 % VC - train ramp for thalamic cells - Christine (for GABAextrasyn or Peptide)
        %% 4: VC - train ramp for thalamic cells
        disp(['Evaluating VC - train tramp for thalamic cells (#4) for ' cellID ' : ' protName ' ...']);
        disp('No code for this analysis has been integrade. Ask Christine for help. If you are Christine, write it now!');
        

    case 11
        %% #11: FI curve (full IV curve or single current injection)'
        % TODO: update save two project data to align with the rest of the VC
        % analysis
        disp(['Evaluating FI curve (#11) for ' cellID ' : ' protName ' ...']);

        resTableName = 'Istep_results';
        if reAnalyze && findProtInProjData(projData,protPath,resTableName)
            disp(['Reanalyzing: ' protName '...']);
        elseif ~reAnalyze && findProtInProjData(projData,protPath,resTableName)
            disp(['Already analyzed, not reanalyzing: ' protName '...'])
        end
        
        %         plotPause = input('Do you want to look at each trace? Yes[1] No[0]...');
        %         if isempty(plotPause)
        %             plotPause = 1;
        %         end
        [spikes_table,Rm,Rs,rheoBase] = calc_FIcurve(protPath,plotPause);
        
        % full FI curve
        if length(unique(spikes_table.I_step))>1 && ~any(isnan(spikes_table.I_step))
            % store results and relavent info in FI_results (struct)
            FI_results.protNum = protNum;
            FI_results.condition = protTable.Condition(ind);
            FI_results.Rm = Rm;
            FI_results.Rs = Rs;
            FI_results.rheoBase = rheoBase;
            FI_results.spikes_table = spikes_table;
            
            % save FI_results (struct) into projData (struct)
            currentResult_str = 'FI_results';
            currentResult = eval(currentResult_str);
            [projData] = saveResult2projData (projData,currentResult_str,currentResult,protNum,cellID,EphysDate,'inputProjData',projData);
            
            % single current injection value
        elseif max(spikes_table.I_step)==min(spikes_table.I_step)
            
            % calculate the mean
            meanArray = nanmean(spikes_table.Variables);
            seArray = nanstd(spikes_table.Variables)./size(spikes_table,1);
            
            meanTable = array2table([meanArray;seArray]);
            meanTable.Properties.VariableNames = spikes_table.Properties.VariableNames;
            meanTable.Properties.RowNames = {'mean','se'};
            
            % store results and relavent info in FI_results (struct)
            Istep_results.protNum = protNum;
            Istep_results.protName = protTable.ProtName(ind);
            Istep_results.condition = protTable.Condition(ind);
            Istep_results.meanSpikesTable = meanTable;
            Istep_results.spikes_table = spikes_table;
            
            % write Istep_results (struct) into projData (struct)
            currentResult_str = 'Istep_results';
            currentResult = eval(currentResult_str);
            [projData] = saveResult2projData (projData,currentResult_str,currentResult,protNum,cellID,EphysDate,'inputProjData',projData);
        end
        
    case 12 % CC - single stim (estim or opto)
        %% 12: CC - single stim
        disp(['Evaluating CC- single stim (#12) for ' cellID ' : ' protName ' ...']);
         
        stimType = protTable.StimType{ind}; %define stim type -- estim, led488, picosprtizer -- should match the ao outputs
        resTableName = 'MGvStim_Results';
        
        % --- find if the data already exists in project data ---
        if (findProtInProjData(projData,protPath,'Results') && reAnalyze ) ... % if we can find the same prot alareay analyzed and that we want to reanalyze
                ||  ~findProtInProjData(projData,protPath,'Results') % or we didn't find it to be alareayd analyzed
            if reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Reanalyzing: ' protName '...']);
            elseif ~reAnalyze && findProtInProjData(projData,protPath,resTableName)
                disp(['Already analyzed, not reanalyzing: ' protName '...'])
            end
            
            % --- load setup ---
            Setup = load(protTable.SetupFile{ind});
%             if plotPause % display and allow changes to setup params
%                 disp(Setup);%display optoSetup to check for any alterations
%                 disp('check for any changes in command window. To edit: optoSetup.xxx = nnn; When finished, write "dbcont"');
%                 keyboard
%             end
            
            % ------------ write polarity ------------
            %check for Recording voltage in protocal table
            [polarity] = writePolarityFromProtTable(protTable,ind);
                        
            % ------- write custom vars 
            extraVars.var1 = protTable.DrugUsed{ind}; extraVars.var1Name = 'DrugUsed';
            extraVars.var2 = protTable.StimLevel(ind); extraVars.var2Name = 'StimLevel';
            extraVars.var3 = protTable.Condition(ind); extraVars.var3Name = 'Condition';
            
            % ----------- analyze each trace and visualize each trace with fit --=
            [resTableName,resTable] = detectStim_CC(cellID,protNum,meta,polarity,Setup,stimType,extraVars,resTableName);
            
            % if user wants to stop and look at all the traces, pause code for user to inspect
            if plotPause; pause;end;   close all; 
            
            % ------------To save struct after analysis into projData
            [projData] = writeProtInProjData(meta, cellID,protNum,resTableName,resTable,'inputProjData',projData);
            
        end
     
    case 13 % CC - hyperpolarizing currents over long period of time
        %% 13: CC - hyperpolarizing currents over long period of time
        disp(['Evaluating CC- hyperpolarizing current currents (#13) for ' cellID ' : ' protName ' ...']);
        disp('No code for this analysis has been integrade. Ask Christine for help. If you are Christine, write it now!');

    case 30
        %% 30: spontaneous data analysis: only calculate Rin, baseline, and peak
        disp(['Evaluating spont data (#30) for ' cellID ' : ' protName ' ...']);
        resTableName = 'spontResults'; 
        if reAnalyze && findProtInProjData(projData,protPath,resTableName)
            disp(['Reanalyzing: ' protName '...']);
        elseif ~reAnalyze && findProtInProjData(projData,protPath,resTableName)
            disp(['Already analyzed, not reanalyzing: ' protName '...'])
        end
            
        
        % calculate the Rin, baseline, and peak for data 
        [resTable] = detectSpontRm(cellID,protNum,meta,resTableName);

         % write into projData (struct)
        [projData] = saveResult2projData (projData,resTableName,resTable,protNum,cellID,EphysDate);
    otherwise
        disp(['No valid analysis found for '  cellID ' : ' protName '. Skipping...']);
 
 %% Extra code to be removed if running properly ----------------------------------
                 %   % --- TO DO: save2projData
           %% [projData] = saveResult2projData (projData,currentResult_str,currentResult,protNum,cellID,EphysDate);

end
end
