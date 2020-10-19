function cluster_rep_table = find_longer_scaffold(Path, myCluster, C, GlobalName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dl=filesep;
    IsolateDir=cellfun(@(c) [Path.Alignment dl 'tree' num2str(myCluster.Num) dl GlobalName  c],C.SeqPlate,'UniformOutput',false);
    AtGFF=cellfun(@(id) GFFAnnotation([id '.gff']), IsolateDir);
    At=arrayfun(@(atg) getData(atg, [1:find(cellfun(@isempty, atg.Attributes), 1, 'first')-1]), AtGFF, 'UniformOutput', false);
    Scaffold=cellfun(@(at) squeeze(split({at.Reference}, '_')), At, 'UniformOutput', false);
    as_node=cellfun(@(as) str2double(as(:,2)), Scaffold, 'UniformOutput', false);
    isoGenes=cellfun(@(at) arrayfun(@(aat) extractBetween(aat, 'Name=',';'), {at.Attributes}, 'UniformOutput', false), At, 'UniformOutput', false);
    for i=1:numel(At)
        isoGenes{i}(cellfun(@isempty, isoGenes{i}))={'group'};
        isoGenes{i}=[isoGenes{i}{:}];
        isoGenes{i}=string(myExtractBefore(isoGenes{i},'_'));
    end
    for i=1:numel(myCluster.Genes)
        cluster_rep_table=[];
        %need to find range of genes in list
        asGenes=string(myExtractBefore(myCluster.Genes{i},'_'))';
        lg=numel(asGenes);
        gn=sum(asGenes~='group');
        %table: seqplate, assembly, cluster number, gene numbe ron isolate, number of genes
        for iso=1:numel(At)
            clusterLoc=find(arrayfun(@(j) sum(isoGenes{iso}(j:j+lg-1)==asGenes ...
                & isoGenes{iso}(j:j+lg-1)~='group')/gn>0.9, 1:numel(isoGenes{iso})-lg));
            for loc=clusterLoc
                if numel(unique(as_node{iso}(loc:loc+lg-1)))==1
                    cluster_rep_table=[cluster_rep_table; iso as_node{iso}(loc) i loc lg];
                end
            end
        end
    end
    if size(cluster_rep_table,1)>1 && numel(unique(cluster_rep_table(:,3)))>1
        cluster_rep_table
        pause
    end


end

