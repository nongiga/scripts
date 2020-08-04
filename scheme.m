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

    
if ~exist('all_alignments20_pangenome.mat', 'file')
    
    pCase=make_pCase(ClusterCase);
    
else
    load('all_alignments20_pangenome.mat', 'pCase'); 
end
end

% % inclusion criteria:
isok = @(isPlas, isTran, isGroup) ~vertcat(isPlas{:}) & ~vertcat(isTran{:}) & ~vertcat(isGroup{:});

delgenes = arrayfun(@(pc, n) [pc.pangenome(horzcat(pc.cGeneID{:}),:) repmat({n}, length(horzcat(pc.cGeneID{:})), 1)], pCase, iClusterCase, 'UniformOutput', 0);
okgenes=cellfun(@(dg) dg(isok(dg(:,3), dg(:,4), dg(:,5)),:), delgenes, 'uniformoutput', 0);
nGenesPerCase=cellfun(@(og) size(og, 1), okgenes);
okgenes=vertcat(okgenes{:});

[ nActual, sigGenes,sgLoc] = calcGeneHist(okgenes(:,1));
figure(1);clf
bar(nActual);
hold on

% randomize
Nrand = 100;
%% calculate number of genes that were deleted multiple times 
nGeneDel_Rand = nan(Nrand,10);
for ir = 1:Nrand
    % Randomize deleted genes:
    randGenes = cell(nC,1);
    for iCase = 1:nC
        pangenome= pCase(iCase).pangenome;
        ok = isok(pangenome(:,3),pangenome(:,4),pangenome(:,5));
        allgenes = pangenome(ok,:); 
        genes = allgenes(:,1);
        p = randperm(numel(genes));
        p = p(1:nGenesPerCase(iCase));
        randGenes{iCase} = genes(p);        
        
    
    end
    randGenes = vertcat(randGenes{:});
    nGeneDel_Rand(ir,:) = calcGeneHist(randGenes);
end

errorbar(1:10,mean(nGeneDel_Rand),std(nGeneDel_Rand));


% is the integral the same?
mean(sum(nGeneDel_Rand*(1:10)',2))
length(okgenes(:,1))

make_figure(okgenes, delgenes, sigGenes, iClusterCase)
%% set the threshold to 4 repetitions and find which genes are those




