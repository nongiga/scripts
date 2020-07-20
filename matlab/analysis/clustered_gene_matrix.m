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
all_descriptions=vertcat(ClusterCase.Description);
all_genes=vertcat(all_genes{:});
all_descriptions=horzcat(all_descriptions{:})';

clean_genes=extractBefore(strcat(all_genes, '_'), '_');
[u, uidx]=unique(clean_genes);


T=table('Size', [numel({ClusterCase.Num}) numel(u)], 'VariableNames', u, ...
    'VariableType', repelem({'uint8'}, numel(u)));
T.Properties.VariableDescriptions=all_descriptions(uidx);
T.Properties.RowNames={ClusterCase.Num};

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];


casenums=arrayfun(@(i) repelem({ClusterCase(i).Num},cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];

is_plasmid=arrayfun(@(i) arrayfun(@(j) repmat(ClusterCase(i).IsPlasmid(j), ClusterCase(i).GeneNum(j),1), 1:numel(ClusterCase(i).GeneNum), 'UniformOutput', false), 1:numel(ClusterCase), 'UniformOutput', false);
is_plasmid=[is_plasmid{:}];
is_plasmid=vertcat(is_plasmid{:});


for i=1:numel(casenums)
    T{casenums{i}, clean_genes{i}}=uint8(is_plasmid(i))+1;
end

save('clustered_genes_table.mat','T');

%% simulate random permutation of numbers in vector
seed=100;
Tmat=table2array(T);
for i=1:height(T)
    seed=seed+1;
    s=rng(seed);
    Tmat(i,:)=Tmat(i, randperm(width(T)));
end

hist(sum(Tmat>0))
