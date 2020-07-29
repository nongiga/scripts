close all; clear all;

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
%geneTable(:,3) = cellfun(@num2str,geneTable(:,3),'UniformOutput',false);
jCase = cat(1,geneTable{:,3});

nGenesPerCase = hist([geneTable{:,3}],1:numel(ClusterCase));

% calculate real number of gene deletions
nGeneDel_Real = genes2hist(geneTable(:,1));

figure(1);clf
bar(nGeneDel_Real);
hold on

% check if the same gene appears twice in a cluster:
isGroup = cellfun(@(s)startsWith(s,'group'),geneTable(:,1));
for iCase = 1:numel(ClusterCase)
    f = jCase == iCase;
    u = unique(geneTable(~isGroup & f,1));
    assert(numel(u)==sum(~isGroup & f));
end
%it does not

Nrand = 100;

%calculate number of genes that were deleted multiple times 
nGeneDel_Rand = nan(Nrand,10);
for ir = 1:Nrand
    % Randomize deleted genes:
    randGenes = cell(nC,1);
    for iCase = 1:nC
        CC = ClusterCase(iCase);
        allgenes = Case(CC.CaseNum).GeneName;
        allgenes = allgenes(:); % coloumn
        p = randperm(numel(allgenes)); 
        p = p(1:nGenesPerCase(iCase));
        randGenes{iCase} = allgenes(p);
    end
    randGenes = cat(1,randGenes{:});
    nGeneDel_Rand(ir,:) = genes2hist(randGenes);
end

errorbar(1:10,mean(nGeneDel_Rand),std(nGeneDel_Rand))


function nGenesWithSpecificNumDel = genes2hist(genes)

% Remove underscore from gene name:
underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), genes(:,1));
isRemove = underscoreLoc>=4 & underscoreLoc<inf;
genes(isRemove,1) = cellfun(@(s) extractBefore(s,'_'),genes(isRemove,1),'UniformOutput',false);

isGroup = cellfun(@(s)startsWith(s,'group'),genes(:,1));

[~,i,j] = unique(genes(~isGroup,1));
nDelPerGene = hist(j,1:numel(i))'; % Number of clusters for each gene
nGenesWithSpecificNumDel = hist(nDelPerGene,1:10);

end