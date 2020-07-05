function Sub=aligning_step(IsolatesNames,pipevar,Path)
%ALIGNING_STEP Summary of this function goes here
%   Detailed explanation goes here
dl=filesep;



save([Path.Output dl 'path'], 'Path')
load([Path.Output dl 'pipevar'], 'pipevar')
load([Path.Output dl 'IsolatesNames','IsolatesNames']);

%load ss_isolates (will at some point have to rewrite this part of the
%program)
%ss_isolates=make_ss_isolates
ss_isolates=readtable(pipevar.same_strain_isolates{:}, 'PreserveVariableNames',true);


%find name of plate
GlobalName=IsolatesNames.RawSubDirName{1};
GlobalName=GlobalName(1:(regexp(GlobalName, '\d','once')-1));

Sub=table([],[],'VariableNames', {'Case','Subdir'});

r=0;

for i=1:height(ss_isolates)
    CaseDir=['tree' num2str(ss_isolates.RandomID(i))];
    CurrentDir=[ Path.Alignment dl CaseDir];
    
    %get all combinations possible of Isolates Plates
    
    SameStrains=erase(ss_isolates{i,[3,6]},'.');
    SameStrains=insertBefore(SameStrains,  cellfun(@(p) find(isletter(p),1), SameStrains), '_');
    
    loc=strcmp(Sub.Case,CaseDir);
    
    if sum(loc)
        Sub.Subdir(loc)={unique([Sub.Subdir{loc}, SameStrains])};
    else
        Sub=[Sub ; {CaseDir, {SameStrains}}];
    end


end

end

