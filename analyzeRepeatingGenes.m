%close all; clear all;

if ~exist('ClusterCase','var')
    load('clusterReports/all_clusters20_report.mat')
end

if ~exist('Case','var')
load('alignmentReports/all_alignments20_report.mat')
end

clc


nC = numel(ClusterCase);
% Build array of all deleted genes:
iClusterCase = 1:nC;
geneTable = arrayfun(@(c,i) [cat(1,c.Genes{:}),cat(1,c.Description{:}),num2cell(i+zeros(numel(cat(1,c.Genes{:})),1))], ClusterCase, iClusterCase, 'UniformOutput',false)';
geneTable = cat(1,geneTable{:});
geneNum=transpose(cat(1,ClusterCase.GeneNum));



%% is it an insertion or deletion?
geneInsert=cat(2,ClusterCase.Insert);
isInsert=arrayfun(@(gi, gn) repmat(gi, gn, 1), geneInsert, geneNum, 'UniformOutput', false);
isInsert=cat(1,isInsert{:}).*1;
isInsert(isInsert==0)=-1;
geneTable(:,4)=num2cell(isInsert);


%% define if it is a plasmid
clusterPlasmid=cat(2,ClusterCase.IsPlasmid);
isPlasmid=arrayfun(@(ci, cn) repmat(ci, cn, 1), clusterPlasmid, geneNum, 'UniformOutput', false);
isPlasmid=cat(1,isPlasmid{:});
geneTable(:,5)=num2cell(isPlasmid);


%% define if it is on a transposon
clusterTra=cat(2,ClusterCase.IsTra);
isTra=arrayfun(@(ci, cn) repmat(ci, cn, 1), clusterTra, geneNum, 'UniformOutput', false);
isTra=cat(1,isTra{:});
geneTable(:,6)=num2cell(isTra);

%% max coverage of gene
maxCov=cat(2,ClusterCase.MaxCov);
maxCov = arrayfun(@(ci,cn) repmat(ci, cn, 1), maxCov, geneNum, 'UniformOutput',false)';
maxCov = cat(1,maxCov{:});
geneTable(:,7)=num2cell(maxCov);


%% 
%geneTable(:,3) = cellfun(@num2str,geneTable(:,3),'UniformOutput',false);
jCase = cat(1,geneTable{:,3});

nGenesPerCase = hist([geneTable{:,3}],1:numel(ClusterCase));

isGroup = cellfun(@(s)startsWith(s,'group'),geneTable(:,1));
%process the gene table to group all undefined similar proteins
geneTableMod=process_gene_table(geneTable, ClusterCase);
% calculate real number of gene deletions
[nGeneDel_Real, sigGenes, sigLoc]= genes2hist(geneTable(~vertcat(geneTable{:,5}) & ~vertcat(geneTable{:,6})));

%get table of significant genes
sigTable=geneTable(startsWith(geneTable(:,1), sigGenes),:);

figure(1);clf
bar(nGeneDel_Real);
hold on


% check if the same gene appears twice in a cluster:
for iCase = 1:numel(ClusterCase)
    f = jCase == iCase;
    u = unique(geneTable(~isGroup & f,1));
    assert(numel(u)==sum(~isGroup & f));
end
%it does not

Nrand = 100;

%% calculate number of genes that were deleted multiple times 
nGeneDel_Rand = nan(Nrand,10);
for ir = 1:Nrand
    % Randomize deleted genes:
    randGenes = cell(nC,1);
    for iCase = 1:nC
        CC = ClusterCase(iCase);
        mC=Case(CC.CaseNum);
        allgenes = mC.GeneName(any(mC.NCov<4,2));
        allgenes = allgenes(:); % coloumn
        p = randperm(numel(allgenes)); 
        p = p(1:nGenesPerCase(iCase));
        randGenes{iCase} = allgenes(p);
    end
    randGenes = cat(1,randGenes{:});
    nGeneDel_Rand(ir,:) = genes2hist(randGenes);
end

errorbar(1:10,mean(nGeneDel_Rand),std(nGeneDel_Rand))



%% create table of function and gene

sigMat=zeros(numel(sigGenes),nC);
for i=1:numel(sigGenes)
    sgLoc=startsWith(sigTable(:,1), sigGenes(i));
    sigMat(i,[sigTable{sgLoc,3}])=[sigTable{sgLoc,4}];
end


% look at the 
figure(2);
IndelSC=[sum(sigMat==1,2) sum(sigMat==-1,2)];
figure(2); h=my_imagesc(IndelSC);
set(gca, 'ytick',[1:numel(sigGenes)],'yticklabel', sigGenes)
set(gca, 'xtick',[1 2],'xticklabel', {'Insertion', 'Deletion'})
colorbar



% dividing the genes to groups
ClusterGenes=arrayfun(@(c) {cat(1,c.Genes{:})}, ClusterCase);
%where are these proteins found? are they on the same strand? or surrounded
%by similar proteins? is it the same phage?
% make a table of cluster, significant genes in it, 


% I don't think we're evaluating this correctly. many of the genes that
% stand out are integration genes and other phage genes. Okay we get it,
% these are phages. But it's not important.

%%

function [nGenesWithSpecificNumDel,  sigGenes, ogSigLoc]= genes2hist(genes)

% Remove underscore from gene name:
correctedGenes=genes(:,1);

underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), correctedGenes);
isRemove = underscoreLoc>=4 & underscoreLoc<inf;
correctedGenes(isRemove) = cellfun(@(s) extractBefore(s,'_'),correctedGenes(isRemove),'UniformOutput',false);

isGroup = cellfun(@(s)startsWith(s,'group'),correctedGenes);

correctedGenes=correctedGenes(~isGroup);

[u,i,j] = unique(correctedGenes,'stable');
nDelPerGene = hist(j,1:numel(i))'; % Number of clusters for each gene
nGenesWithSpecificNumDel = hist(nDelPerGene,1:10);

%get all the locations in genes of significant genes
sigLoc=find(hist(j, unique(j))>3);
sigGenes=u(sigLoc);

gLoc=find(~isGroup);
ogSigLoc=gLoc(i(sigLoc));

end