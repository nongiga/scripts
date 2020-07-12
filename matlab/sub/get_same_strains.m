function [SameStrains, GlobalName]=get_same_strains(ss_isolates_file, IsolatesNames, GlobalName, j)

ss_isolates=readtable(ss_isolates_file, 'PreserveVariableNames',true);
dl=filesep;

%find name of plate


SameStrains=table([],[],'VariableNames', {'Case','seqPlates'});

for i=1:height(ss_isolates)
    RID=num2str(ss_isolates.RandomID(i));
    
    %get all combinations possible of Isolates Plates
    
    seqPlates=erase(ss_isolates{i,[3,6]},'.');
    seqPlates=insertBefore(seqPlates,  cellfun(@(p) find(isletter(p),1), seqPlates), '_');
    
    %take out same strains not in isolates list
    seqPlates=seqPlates(arrayfun(@(s) any(endsWith(IsolatesNames.RawSubDirName, strcat(GlobalName, s))), seqPlates));
    
    loc=strcmp(SameStrains.Case,RID);
    if sum(loc)
        SameStrains.seqPlates(loc)={unique([SameStrains.seqPlates{loc}, seqPlates])};
    else
        SameStrains=[SameStrains ; {RID, {seqPlates}}];
    end


end
SameStrains=SameStrains(cellfun(@length, SameStrains.seqPlates)>1,:);

InstructionsRef=repelem(j,height(SameStrains))';
SameStrains=addvars(SameStrains, InstructionsRef);





