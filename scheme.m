%% Data base sketch

if ~exist('ClusterCase','var')
    load('clusterReports/all_clusters20_report.mat')
end




nC = numel(ClusterCase);
iClusterCase = 1:nC;
ClusterNum=[ClusterCase.CaseNum];
% for each case iCase
%%make array of cell tables name pangenome
if ~exist('pCase', 'var')
    if ~exist('Case','var')
        load('alignmentReports/all_alignments20_report.mat')
    end
    
if ~exist('all_alignments20_pangenome.mat', 'file')
    %[Case.pangneome]=NaN; [Case.cGeneIDs]=NaN;
    pCase=struct('pangenome', {}, 'cGeneID',{});
    for iCase=iClusterCase
        % pangenome = {'gene', 'description', isPlasmid, isTranspsone, isGroup, hasName} % table of ~4000 genes
        mC=Case(ClusterNum(iCase));

        isPlasmid= num2cell(any(contains(lower(mC.GeneDescription), 'plasmid') | mC.NCov>4, 2));
        isTransposon=num2cell(any(contains(lower(mC.GeneDescription), 'transpos') | startsWith(mC.GeneName, 'tra'), 2));
        isGroup = num2cell(arrayfun(@(s)startsWith(s,'group'),mC.GeneName));

        givenName=give_gene_name(mC.GeneName, mC.GeneDescription);
        hasName=num2cell(~arrayfun(@(s)startsWith(s,'group'),givenName));

        pCase(iCase).pangenome=[ mC.GeneName(:), mC.GeneDescription, isPlasmid, isTransposon, isGroup,hasName, givenName];

        %if this case contains deletions/insertions
        genelocs=ClusterCase(iCase).Loc;
        pCase(iCase).cGeneID=arrayfun(@(i) [genelocs(i,1):genelocs(i,2)],1:size(genelocs,1), 'UniformOutput',0);

    end
    disp('saving')
    save('all_alignments20_pangenome.mat', 'pCase','-v7.3')
    clear Case
    disp('finished saving')
else
    load('all_alignments20_pangenome.mat', 'pCase'); 
end
end

% % inclusion criteria:
isok = @(isPlas, isTran, isGroup) ~vertcat(isPlas{:}) & ~vertcat(isTran{:}) & ~vertcat(isGroup{:});

delgenes = arrayfun(@(pc) pc.pangenome(horzcat(pc.cGeneID{:}),:) , pCase, 'UniformOutput', 0);
okgenes=cellfun(@(dg) dg(isok(dg(:,3), dg(:,4), dg(:,5)),:), delgenes, 'uniformoutput', 0);
nGenesPerCase=cellfun(@(og) size(og, 1), okgenes);
okgenes=vertcat(okgenes{:});
%vertcat(delgenes{:}
%??? Case, ClterGene, isok)

nActual = calcGeneHist(okgenes(:,1));
figure(1);clf
bar(nActual);
hold on

% randomize
Nrand = 100;
%% calculate number of genes that were deleted multiple times 
nGeneDel_Rand = nan(Nrand,10);
for ir = 1:Nrand
    % Randomize deleted genes:
    randGenes = cell(nC,1); ClusterNum
    for iCase = 1:nC
        pangenome= pCase(iCase).pangenome;
        ok = isok(pangenome(:,3),pangenome(:,4),pangenome(:,5));
        allgenes = pangenome(ok,:); 
        genes = allgenes(:,1);
        p = randperm(numel(genes));
        p = p(1:nGenesPerCase(iCase));
        randGenes{iCase} = allgenes(p);        
        
    
    end
    randGenes = [randGenes{:}];
    nGeneDel_Rand(ir,:) = calcGeneHist(randGenes);
end

errorbar(1:10,mean(nGeneDel_Rand),std(nGeneDel_Rand))






function NumGenesDeletedNtimes = calcGeneHist(geneNames)
% geneNames: names of genes deleted in each Case
% {{'gene1','gene2',...}, {'gene1','gene5'}}

% Remove underscore from gene name:
    underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), geneNames);
    isRemove = underscoreLoc>=4 & underscoreLoc<inf;
    geneNames(isRemove) = cellfun(@(s) extractBefore(s,'_'),geneNames(isRemove),'UniformOutput',false);

    [u,i,j] = unique(geneNames,'stable');
    nDelPerGene = hist(j,1:numel(i))'; % Number of clusters for each gene
    NumGenesDeletedNtimes = hist(nDelPerGene,1:10);
    

end



function givenNames=give_gene_name(geneName, geneDesc)

% give names based on NAME in descriptions
    isGroup = cellfun(@(s)startsWith(s,'group'),geneName);
    isHypot=contains(geneDesc, 'hypothetical protein');
    nNameyDesc=find(isGroup & ~isHypot);
    lastWord=cellfun(@extractLast, geneDesc(nNameyDesc));
    LWlen=cellfun(@length, lastWord);
    isGeneName=LWlen > 1 & ...
        ((cellfun(@(s) any(regexp(s ,'[0-9]')), lastWord)) | ...
        (LWlen >= 3 & LWlen <=4)) & ...
        ~cellfun(@(s) contains(s, 'ase'), lastWord);
    geneName(nNameyDesc(isGeneName))=lastWord(isGeneName);
    
    
    % clean name underscores and whatnot
    isGroup = cellfun(@(s)startsWith(s,'group'),geneName);
    underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), geneName);
    isRemove = underscoreLoc>=4 & underscoreLoc<inf & ~isGroup;
    geneName(isRemove) = cellfun(@(s) extractBefore(s,'_'),geneName(isRemove),'UniformOutput',false);

    %% if same description as a different, named protein, apply it to the unnamed proteins
    [u_desc, ~, ud_idx]=unique(geneDesc);
    for i=1:max(ud_idx)
        sameDesc = ud_idx==i;
        if nnz(sameDesc) > nnz(sameDesc & isGroup)
            newGeneName=unique(lower(geneName(sameDesc & ~isGroup)));
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
                    geneName(sameDesc & isGroup)= clearNewGeneName(1);
                end

            else
                geneName(sameDesc & isGroup)= newGeneName;
            end

        end
    end

    givenNames=geneName;
    
    function s=extractLast(string)
    sp=split(string, [" ", "-"]);
    s=sp(end);
end
end
