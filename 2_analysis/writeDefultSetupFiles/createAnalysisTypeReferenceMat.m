%% this code defines the analysis type and setup files

% Created by CJL, 11/28/2022

% KEY FOR ANALYSIS INDEX
% VC analyis:
% 1: minimal e stimulation analysis (detectSinglePeak.m) - no charge calc
% 2: single opto stimulation analysis (detectOptoStim.m)
% 3: train eStim 
% 4: train opto stim
% 5: train ramp (Christine's thalamic cells)

% CC analysis 
% 11: full FI curve (full IV curve or single current injection)
% 12: hyperpolarizing currents over long period of time 
% 13: opto train in CC



%% Create an empty table with 3 columns: each of the following sections can add to this overall analysis type table 

emptyTable = table('Size',[0,7],...%create an empty table with 3 columns
        'VariableTypes',{'double','string','cell','cell','cell','cell','cell'},... % specify variable type
        'VariableNames',["analysisIndex" "analysisName" "SetupFile" "otherParams" "paramNotes" "checkStr" "excludeStr"]);  % name the headers of these variables   
analysisType = emptyTable;

%% 1. VC minimal evoked
Min = emptyTable;
Min.analysisIndex(1) = 1; 
Min.analysisName(1) = "min evoked";

Min.SetupFile(1) = {'\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\monoSetup_default.mat'};
Min.otherParams (1) = {{'RecordingVoltage','StimType','StimLevel'}};
Min.paramNotes (1) = {{'must be number','must match the ao list: ex - estim, led488','number or str'}};

Min.checkStr(1) = {{'min'}};
Min.excludeStr(1) = {{'DNU','find','abort'}};
thisTable = Min; %so that we don't need to keep changing the saving and updating part.

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable
%% 2. VC single evoked -- could be either Opto or estim
Single = emptyTable;
Single.analysisIndex(1) = [2]; 
Single.analysisName(1) = ["single evoked VC"];
Single.SetupFile(1) = {{'\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\optoSetup_default.mat',...
    '\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\eStimSetup_default.mat'}};
Single.otherParams (1) = {{'RecordingVoltage','Polarity','Condition','StimType','StimLevel','DrugUsed'}};
Single.paramNotes (1) = {{'must be number','number: 1 or -1','user defined string','string; must match the ao list: ex - estim, led488','number or str','str'}};

Single.checkStr(1) = {{'stim','x threshold','mA'}};
Single.excludeStr(1) = {{'apply','DNU','find','abort'}};
thisTable = Single;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable

%% 3. VC train evoked (n>1 number of stims)  -- could be either Opto or estim

Train = emptyTable;
Train.analysisIndex(1) = [3]; 
Train.analysisName(1) = ["train evoked"];
Train.SetupFile(1) = {'\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\trainSetup_default.mat'};
Train.otherParams (1) = {{'RecordingVoltage','Polarity','TrainPulse','Condition','StimType','StimLevel','DrugUsed'}};
Train.paramNotes (1) = {{'must be number','number: 1 or -1','number in Hz','user defined string','string; must match the ao list: ex - estim, led488','number or str','str'}};
Train.checkStr(1) = {{'train','pairedPulse','Hz'}};
Train.excludeStr(1) = {{'DNU','abort'}};
thisTable = Train;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable


%% 4. Opto train ramp - thalamic
thalamicOpto = emptyTable;
thalamicOpto.analysisIndex(1) = [4]; 
thalamicOpto.analysisName(1) = "Opto train ramp - thalamic";

thalamicOpto.SetupFile(1) = {''};
thalamicOpto.otherParams (1) = {{'FreqRamp','DrugUsed','Condition','RecordingVoltage'}};
thalamicOpto.paramNotes (1) = {{'must be number','str','user defined string','must be number'}};

thalamicOpto.checkStr(1) = {{'Opto train ramp'}};
thalamicOpto.excludeStr(1) = {{'MGv','abort'}};

thisTable = thalamicOpto;
% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable


%% 11. CC - FI curve (full IV curve or single current injection)
FIcurve = emptyTable;
FIcurve.analysisIndex(1) = [11]; 
FIcurve.analysisName(1) = ["FI curve"];

FIcurve.SetupFile(1) = {''};
FIcurve.otherParams (1) = {{'Condition','DrugUsed'}};
FIcurve.paramNotes (1) = {{'user defined str','str'}};
FIcurve.checkStr(1) = {{'IV curve','FI curve'}};
FIcurve.excludeStr(1) = {{'check','find','abort','DNU'}};

thisTable = FIcurve;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable

%% 12. CC - single stim 
SingleCC = emptyTable;
SingleCC.analysisIndex(1) = [12]; 
SingleCC.analysisName(1) = ["single evoked CC"];
SingleCC.SetupFile(1) = {{'\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\optoSetup_CC_default.mat',...
    '\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\eStimSetup_CC_default.mat'}};
SingleCC.otherParams (1) = {{'HoldCurent','Polarity','Condition','StimType','StimLevel','DrugUsed'}};
SingleCC.paramNotes (1) = {{'must be number','number: 1 or -1','number: 1 or -1','user defined string','string; must match the ao list: ex - estim, led488','number or str','str'}};

SingleCC.checkStr(1) = {{'stim','mA'}};
SingleCC.excludeStr(1) = {{'apply','DNU','find'}};
thisTable = SingleCC;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable

%% 13. CC - train stim 

%% 14. CC - Vrest w Hyperpolarizing Current Steps (Christine peptide)
thalamicPeptide = emptyTable;
thalamicPeptide.analysisIndex(1) = [13]; 
thalamicPeptide.analysisName(1) = ["Vrest w Hyperpolarizing Current Steps"];

thalamicPeptide.SetupFile(1) = {''};
thalamicPeptide.otherParams (1) = {{'Condition','DrugUsed'}};
thalamicPeptide.paramNotes (1) = {{'user defined str','str'}};

thalamicPeptide.checkStr(1) = {{'current steps'}};
thalamicPeptide.excludeStr(1) = {{'MGv'}};

thisTable = thalamicPeptide;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable


%% 30. CC spont data
spont = emptyTable;
spont.analysisIndex(1) = [11]; 
spont.analysisName(1) = ["Spont Rin"];

spont.SetupFile(1) = "";
spont.otherParams (1) = {{'Condition','DrugUsed'}};
spont.checkStr(1) = {'spont'};
spont.excludeStr(1) = {'DNU'};

thisTable = spont;

% ----------- save thisTblae into the overall analysisType table 
analysisType = [analysisType;thisTable];
clear thisTable


%% save the analysisType variable into the common ephys data analysis folder on apollo
analysisType.Properties.Description = ['Please check the analysisType table loaded in the variable window and make approperiate changes to checkStr column. \n'...
'Each strings in the checkStr column will be compared to the names of every protocal folder; \n'...
'if the protocal folder name contains the checkStr for a given analysis type, then protocal belongs to this type of analysis.'];
fileName = fullfile("\\apollo\research\ENT\Takesian Lab\Ephys Data Analysis\0_settings\","analysisTypeReference.mat");
save (fileName,'analysisType');
disp('new analysisType .mat structure has been saved onto Apollo (Ephys Data Analysis\0_settings');