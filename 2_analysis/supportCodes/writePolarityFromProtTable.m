function [polarity] = writePolarityFromProtTable(protTable,ind)
%check for polarity in protocal table given the ind (row of the protTable)

%%
if  any(strcmp('Polarity', protTable.Properties.VariableNames)) % if protocal table has the column "RecordingVoltage"
    if ~isempty(protTable.Polarity(ind)) %check if this entry is a number
        polarity = protTable.Polarity(ind);
        % if not in protocal table or data entry not numeric, ask for manual input
    else
        polarity = input('EPSCs = -1; IPSCs = 1...');
    end
elseif any(strcmp('RecordingVoltage', protTable.Properties.VariableNames)) % if protocal table has the column "RecordingVoltage"
    if isnumeric(protTable.RecordingVoltage(ind)) %check if this entry is a number
        if protTable.RecordingVoltage(ind) > 0
            polarity = 1;
        else
            polarity = -1;
        end
    else % if not in protocal table or data entry not numeric, ask for manual input
        polarity = input('EPSCs = -1; IPSCs = 1...');
    end
    
end