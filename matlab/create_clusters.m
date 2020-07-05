
close all
clear all
dl=filesep;
load('gap_data');

disp("create clusters")
for ins=1:h
    ClusterCase=[];
    m=moptions{pipevar.report_multi(ins)+1};
    id=[ num2str(pipevar.bp(ins)) m  ];   
    if ~pipevar.recreate_reports(ins) && exist([Path.Clusters  dl 'all_alignments' id '_cluster_report.mat'],'file')
        continue
    end
    
    load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')        

    
    for iCase = 1: numel(Case)

        C=Case(iCase);        
        [mn, mnidx]=min(C.NCov, [],2);
        [mx, mxidx]=max(C.NCov, [], 2);
        del=mx>0.5 & mn<0.05;
        sum(del)
        if sum(del)
            load([Path.Reports dl 'tree' num2str(C.Num) '_' id] )

            %define the beginning & end of clusters (inclusive)
            cluster_loc=find(arrayfun(@(i) del(i)~=del(i-1), 2:numel(del)));
            if mod(size(cluster_loc, 2),2)~=0
                cluster_loc=[cluster_loc length(del)];
            end
            myCluster=struct;

            myCluster.Loc=reshape(cluster_loc, 2, [])';
            myCluster.Loc(:,1)= myCluster.Loc(:,1)+1;

            myCluser.Num=C.Num;
            myCluster.GeneNum=myCluster.Loc(:,2)-myCluster.Loc(:,1)+1;
            Locs=arrayfun(@(i) myCluster.Loc(i,1):myCluster.Loc(i,2), 1:length(myCluster.GeneNum), 'UniformOutput', 0);
            myCluster.Length=arrayfun(@(i) sum(myCase.GeneLength(Locs{i})), 1:length(Locs));
            myCluster.Genes=arrayfun(@(i) C.GeneName(Locs{i}), 1:length(Locs), 'UniformOutput',false);
            myCluster.Insert=[C.Date{mxidx(myCluster.Loc(:,1))}]>[C.Date{mnidx(myCluster.Loc(:,1))}];
            myCluster.MaxCov=arrayfun(@(i) mean(mx(Locs{i})), 1:length(Locs));
            myCluster.MinCov=arrayfun(@(i) mean(mn(Locs{i})), 1:length(Locs));
            myCluster.IsPlasmid=myCluster.MaxCov>10;

            save([Path.Clusters dl 'tree' num2str(C.Num) '_' id], 'myCluster' )
            ClusterCase=[ClusterCase myCluster];
        end
        
    end
    save([Path.Clusters  dl 'all_alignments' id '_cluster_report.mat'], 'ClusterCase');
end






