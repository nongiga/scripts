function InstructionLines = get_instructions

%get instructions from excel file, turn them into structure


%read and clean table into structure
t= readtable('instructions.xlsx') ;
t.Properties.VariableNames=strrep(t.Properties.VariableNames, "x", "");

%ensure there is an isolates list file
fnm='isolates_list.xlsx';
if ~ismember('isolates_file_name',t.Properties.VariableNames)
    t.isolates_file_name = repelem(fnm, height(t),1);
end

InstructionLines=t(t.Used==1,:);




return