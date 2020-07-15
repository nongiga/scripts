
close all
clear all
dl=filesep;
load('gap_data');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';
disp("create clusters")
h=1;
for ins=1:h
    ClusterCase=[];
    m=moptions{pipevar.report_multi(ins)+1};
    id=[ num2str(pipevar.bp(ins)) m  ];   

    
    load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')        

    
    for iCase = 1: numel(Case)

        
        C=Case(iCase);        
        [mn, mnidx]=min(C.NCov, [],2);
        [mx, mxidx]=max(C.NCov, [], 2);
        del=mx>0.5 & mn<0.05;
        sum(del)
        if sum(del) && numel(del)<6000 && exist([Path.Reports dl 'tree' num2str(C.Num) '_' id  '.mat'], 'file')
            load([Path.Reports dl 'tree' num2str(C.Num) '_' id  '.mat'] )
            
            
            gs=cellfun(@(g) g(1:find(g=='_', 1,'last')-1), myCase.GeneSource, 'UniformOutput', false);
            %define the beginning & end of clusters (inclusive)
            cl=get_cluster_location(mn,mx,del, myCase.AssemblyNum, gs);
            
            if isempty(cl),     continue; end
            
            myCluster=struct;
            myCluster.Loc= cl;
            myCluster.CaseNum=iCase;
            myCluster.Num=C.Num;
            myCluster.GeneNum=myCluster.Loc(:,2)-myCluster.Loc(:,1)+1;
            Locs=arrayfun(@(i) myCluster.Loc(i,1):myCluster.Loc(i,2), 1:length(myCluster.GeneNum), 'UniformOutput', 0);
            myCluster.Length=arrayfun(@(i) sum(myCase.GeneLength(Locs{i})), 1:length(Locs));
            myCluster.Genes=arrayfun(@(i) C.GeneName(Locs{i}), 1:length(Locs), 'UniformOutput',false);
            myCluster.Description=arrayfun(@(i) myCase.GeneDescription(Locs{i}), 1:length(Locs), 'UniformOutput',false)';
            myCluster.Assembly=arrayfun(@(i) myCase.AssemblyNum(myCluster.Loc(i,1)), 1:length(Locs));
            myCluster.Insert=[C.Date{mxidx(myCluster.Loc(:,1))}]>[C.Date{mnidx(myCluster.Loc(:,1))}];
            myCluster.MaxCov=arrayfun(@(i) mean(mx(Locs{i})), 1:length(Locs));
            myCluster.MinCov=arrayfun(@(i) mean(mn(Locs{i})), 1:length(Locs));
            myCluster.IsPlasmid=myCluster.MaxCov>10;
            isTra=arrayfun(@(i) any(contains(myCluster.Genes{i}, 'tra')) | ...
                any(contains(myCluster.Description{i}, 'transposon')), 1:length(myCluster.Genes), 'UniformOutput', 0);
            myCluster.isTra=[isTra{:}];
            save([Path.Clusters dl 'tree' num2str(C.Num) '_' id], 'myCluster' )
            ClusterCase=[ClusterCase myCluster];
        end
        
    end
    save([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase');
end






