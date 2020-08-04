function [ NumGenesDeletedNtimes, sigGenes, ogSigLoc] = calcGeneHist(geneNames)
% geneNames: names of genes deleted in each Case
% {{'gene1','gene2',...}, {'gene1','gene5'}}

% Remove underscore from gene name:
    underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), geneNames);
    isRemove = underscoreLoc>=4 & underscoreLoc<inf;
    geneNames(isRemove) = cellfun(@(s) extractBefore(s,'_'),geneNames(isRemove),'UniformOutput',false);

    [u,i,j] = unique(geneNames,'stable');
    nDelPerGene = hist(j,1:numel(i))'; % Number of clusters for each gene
    NumGenesDeletedNtimes = hist(nDelPerGene,1:10);
    
    
    sigLoc=find(hist(j, unique(j))>3);
    sigGenes=u(sigLoc);
    fi=find(i);
    ogSigLoc=fi(sigLoc);
    
end