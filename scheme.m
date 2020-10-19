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

 % Name(1),Description(2), isPlasmid(3), isTransposon(4), isGroup(5)
% ,hasName(6), givenName(7), isPlasmidCorrect(8), isGenomic(9),isPhage(10),
 % casenum(11), isInsert(12), clusterNum(13)

% % inclusion criteria:
isok = @(isPlas, isTran, isGroup, isGenomic) ~vertcat(isPlas{:}) & ~vertcat(isTran{:}) & ~vertcat(isGroup{:}) & vertcat(isGenomic{:});

delgenes = arrayfun(@(pc, n) [pc.pangenome(horzcat(pc.cGeneID{:}),:) repmat({n}, length(horzcat(pc.cGeneID{:})), 1) num2cell(pc.cInsert') num2cell(pc.cNum') num2cell(pc.cPhage')], pCase, iClusterCase, 'UniformOutput', 0);

%this line includes tra and does not limit scaffold length
%okgenes=cellfun(@(dg) dg(isok(dg(:,8), dg(:,4), dg(:,5), repmat({1}, size(dg,1),1)),:), delgenes, 'uniformoutput', 0);
okgenes=cellfun(@(dg) dg(isok(dg(:,8), dg(:,4), dg(:,5), dg(:,9)),:), delgenes, 'uniformoutput', 0);
nGenesPerCase=cellfun(@(og) size(og, 1), okgenes);
okgenes=vertcat(okgenes{:});
delgenes=vertcat(delgenes{:});
[ nActual, sigGenes,sgLoc] = calcGeneHist(okgenes(:,1), 2);
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
        ok = isok(pangenome(:,8),pangenome(:,4),pangenome(:,5), pangenome(:,9));
        allgenes = pangenome(ok,:); 
        genes = allgenes(:,1);
        p = randperm(numel(genes));
        p = p(1:nGenesPerCase(iCase));
        randGenes{iCase} = genes(p);        
        
    
    end
    randGenes = vertcat(randGenes{:});
    nGeneDel_Rand(ir,:) = calcGeneHist(randGenes,2);
end

errorbar(1:10,mean(nGeneDel_Rand),std(nGeneDel_Rand));


% is the integral the same?
mean(sum(nGeneDel_Rand*(1:10)',2))
length(okgenes(:,1))

%sigGenes=0

make_figure(okgenes, delgenes, sigGenes, iClusterCase)

