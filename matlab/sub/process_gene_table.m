function geneTable=process_gene_table(geneTable, ClusterCase)
%% set up variables
Path.Reports='alignmentReports';
Path.Blast='blast';
dl=filesep; id='20';

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];

casenums=arrayfun(@(i) repelem(i,cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];
all_locs=[]; all_sequences=[];

entry_name=join( [geneTable(:,1), cellstr(num2str(casenums'))], '-') ;

%% process geneTable to give more group names and descriptions

% give names based on NAME in descriptions
isGroup = cellfun(@(s)startsWith(s,'group'),geneTable(:,1));
isHypot=contains(geneTable(:,2), 'hypothetical protein');
nNameyDesc=find(isGroup & ~isHypot);
lastWord=cellfun(@extractLast, geneTable(nNameyDesc,2));
LWlen=cellfun(@length, lastWord);
isGeneName=LWlen > 1 & ...
    ((cellfun(@(s) any(regexp(s ,'[0-9]')), lastWord)) | ...
    (LWlen >= 3 & LWlen <=4)) & ...
    ~cellfun(@(s) contains(s, 'ase'), lastWord);
geneTable(nNameyDesc(isGeneName), 1)=lastWord(isGeneName);

% give description based on name
nameSearch=geneTable(~isGroup & isHypot,1);
isRemove=contains(nameSearch,'_');
nameSearch(isRemove)=extractBefore(nameSearch(isRemove), '_');
descSearch=cellfun(@(ns) find(startsWith(geneTable(:,1),ns) & ~isHypot) , nameSearch,'UniformOutput', false);
descLocs=cellfun(@(ds) ds(1), descSearch);
geneTable(~isGroup & isHypot,2)=geneTable(descLocs, 2);

%% remove underscore/number
isGroup = cellfun(@(s)startsWith(s,'group'),geneTable(:,1));

underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), geneTable(:,1));
isRemove = underscoreLoc>=4 & underscoreLoc<inf & ~isGroup;
geneTable(isRemove,1) = cellfun(@(s) extractBefore(s,'_'),geneTable(isRemove,1),'UniformOutput',false);

%% if same description as a different, named protein, apply it to the unnamed proteins

%loop through unique descriptions
% 4_3: same name, different description

[u_desc, ~, ud_idx]=unique(geneTable(:,2));
for i=1:max(ud_idx)
    sameDesc = ud_idx==i;
    if nnz(sameDesc) > nnz(sameDesc & isGroup)
        newGeneName=unique(lower(geneTable(sameDesc & ~isGroup,1)));
        %issue: the names cannot be applied because there is a number of
        %proteins with same description, different names
        if length(newGeneName)>1
            %if all these names are the same if you remove anything
            %following an underscore...
            isRem=contains(newGeneName,'_');
            clearNewGeneName=newGeneName;
            clearNewGeneName(isRem)=extractBefore(clearNewGeneName(isRem), '_');
            clearNewGeneName=erase(clearNewGeneName, ["orf", "gp"]);
            clearNewGeneName(ismember(clearNewGeneName, 'np'))={'n'};
            
            sameOneLetter=length(unique(clearNewGeneName))==1;
            %or start in the same 3 letters
            if all(cellfun(@length, newGeneName)>=3)
                firstThree=cellfun(@(g) g(1:3), newGeneName, 'UniformOutput', 0);
                startSame=length(unique(firstThree))==1;
            else
                startSame=0;
            end
            
            if sameOneLetter || startSame || length(unique(clearNewGeneName))==1
                geneTable(sameDesc & isGroup,1)= clearNewGeneName(1);
            else
                disp(newGeneName)
            end
            
        else
            geneTable(sameDesc & isGroup,1)= newGeneName;
        end
        
    end
end

%% make sure proteins with certain name and description do not overlap 

%% get sequence of all genes
blastfaa=[Path.Blast dl 'deleted_sequences.fa'];
if ~exist(blastfaa, 'file')
    for iCC=1:numel(ClusterCase)
        CC=ClusterCase(iCC);
        load([Path.Reports dl 'tree' CC.Num '_' id '.mat'])

        locs = arrayfun(@(i) CC.Loc(i,1):CC.Loc(i,2), ...
            1:size(CC.Loc, 1), 'UniformOutput', 0);
        locs=[locs{:}];
        all_locs= [all_locs locs];
        cluster_sequences=myCase.Sequence(locs);
        all_sequences=[ all_sequences; cluster_sequences];

    end

    % make fasta files and run cd-hit
    
    mkdir(Path.Blast);
    p_seq=nt2aa(all_sequences);
    blastfaa=[Path.Blast dl 'deleted_sequences.fa'];
    fastawrite(blastfaa, entry_name, p_seq);
end
if ~exist('blast/deleted_sequences.fa.groups', 'file')
    
    system('iterative_cdhit -m blast/deleted_sequences.fa -f blast/deleted_sequences.fa -c blast/deleted_sequences -v');
end

%% group all proteins based on cluster

ClusterT=readtable('blast/deleted_sequences.fa.groups', 'FileType', 'text');
%iterate through each cluster
for i=1:height(ClusterT)
    gene_loc=ismember(entry_name, ClusterT{i,:});
    c= ismember(geneTable(gene_loc,2), {'hypothetical protein', 'putative protein'});
    % if all the proteins are hypothetical change desc
    if all(c)
        geneTable(gene_loc,[1 2])={['cluster' num2str(i)]};
    %if one or more of the genes is defined but the rest are not
    elseif sum(c)
        % apply the description to all hypothetical proteins
        %divide descriptions to 
        descriptions=geneTable(gene_loc & ~isHypot,2);
        geneTable(gene_loc & isHypot)=descriptions(1);
        
    end
    
end


function s=extractLast(string)
    sp=split(string, [" ", "-"]);
    s=sp(end);
end

end