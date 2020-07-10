
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
% 
% for iCase=1:numel(ClusterCase)
%     
%     Cl=ClusterCase(iCase);
%     load([Path.Reports dl 'tree' num2str(C.Num) '_20.mat'],'myCase' )
%     C=myCase;
% 
%     del=arrayfun(@(i) Cl.Loc(i,1):Cl.Loc(i,2), 1:size(Cl.Loc,1), 'UniformOutput',false);
%     del=[del{:}];
%     
%     genes=vertcat(Cl.Genes{:});
%     ass_num=cellstr(num2str(myCase.AssemblyNum(del)));
%     genes=join([genes, ass_num]);
%     
%     
%     [mn, mnidx]=min(C.NCov, [],2);
%     [mx, mxidx]=max(C.NCov, [], 2);
%     
%     plot(myCase.GeneLength, 
%     
% end


plot_clusters(Path, moptions, pipevar, m, ClusterCase, Case)



