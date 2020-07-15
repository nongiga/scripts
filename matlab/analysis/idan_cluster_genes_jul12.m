clear all; close all;
%% set-up
dl=filesep;
load('gap_data', 'moptions', 'pipevar', 'Path');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

ins=1;
m=moptions{pipevar.report_multi(ins)+1};
id=[ num2str(pipevar.bp(ins)) m  ];  

load([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase')
load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')

%% get all gene names
all_genes=[ClusterCase.Genes];
all_genes=vertcat(all_genes{:});
%how many genes have the name "group"? (meaning they're undefined
undefined=contains(all_genes, 'group_');
sum(undefined)/numel(all_genes)

%% get the sequences of these genes
all_sequences={};
all_genes={};
for iCase=1:numel(ClusterCase)
        
    Cl=ClusterCase(iCase);
    load([Path.Reports dl 'tree' Cl.Num '_' id '.mat'], 'myCase');
    for i=1:size(Cl.Loc,1)
        all_sequences=[all_sequences; myCase.Sequence(Cl.Loc(i,1):Cl.Loc(i,2))];
        all_genes=[all_genes; cellfun(@(c)[c '_' Cl.Num],vertcat(Cl.Genes{:}),'uni',false)];
    end
end

%get sequences for unidentified groups
all_sequences(undefined);

fastawrite('undefined.fasta', all_genes(undefined), all_sequences(undefined))

