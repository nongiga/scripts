clear all; close all;
dl=filesep;
load('gap_data', 'moptions', 'pipevar', 'Path');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

ins=1;
m=moptions{pipevar.report_multi(ins)+1};
id=[ num2str(pipevar.bp(ins)) m  ];  
load([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase')
% 
% 
load([Path.Reports  dl 'all_alignments' id 'trun_report.mat'], 'Case')      
% 
% Inserts=arrayfun(@(c) mean(c.Insert), ClusterCase);
% sum(Inserts==1)
% sum(Inserts==0)
% sum(Inserts>0 & Inserts<1)

% I am seeing (anecdotally) a lot of short genes with very low coverage. what is the difference between
% min(NCOV) and max(NCOV)?



plot_clusters(Path,pipevar, m, ClusterCase, Case, ins)



