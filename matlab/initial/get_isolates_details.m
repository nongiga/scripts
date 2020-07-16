function IsolatesNames = get_isolates_details(fnm,i)


t = readtable(fnm, 'ReadVariableNames', true) ;
t.Properties.VariableNames=strrep(t.Properties.VariableNames, "x", "");
t=(t(t.Used==1,:));
InstructionsRef=repelem(i,height(t))';
IsolatesNames=addvars(t, InstructionsRef);


return