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

[u, uidx]=unique(extractBefore(strcat(all_genes, '_'), '_'));

T=table();
T.Properties.VariableNames=u;
T.Properties.VariableDescriptions=all_descriptions(uidx);
T.Properties.RowNames={ClusterCase.Num};