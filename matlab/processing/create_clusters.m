
% close all
% clear all
dl=filesep;
load('gap_data');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';
disp("create clusters")
h=1;

phage_genes={'bet','Alpha','S_','N_','C_','4_','3_','gin','cox',...
    'B_','Q_','O_','L_','X_','Y_','R_','V_','repC','bor','ant','Alpha',...
    'ninB','cre','pol', 'cI','C2','rep', 'MOD',...
    'P_','cII','gam','bet','exo', 'yokD', 'gag'};

plasmid_genes={'mbe', 'rop'};

for ins=1:h
    ClusterCase=[];
    m=moptions{pipevar.report_multi(ins)+1};
    id=[ num2str(pipevar.bp(ins)) m  ];   

%     load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')        

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
            
            myCluster.Genes=arrayfun(@(i) C.GeneName(Locs{i}), 1:length(Locs), 'UniformOutput',false);
            myCluster.Description=arrayfun(@(i) myCase.GeneDescription(Locs{i}), 1:length(Locs), 'UniformOutput',false)';
            
            myCluster.Source=arrayfun(@(i) myCase.GeneSource(myCluster.Loc(i,1)), 1:length(Locs));
            
            myCluster.Assembly=arrayfun(@(i) myCase.AssemblyNum(myCluster.Loc(i,1)), 1:length(Locs));
            
            myCluster.AssemblyStart=arrayfun(@(i) myCase.GeneStart(myCluster.Loc(i,1)), 1:length(Locs));
            myCluster.AssemblyEnd=arrayfun(@(i) myCase.GeneEnd(myCluster.Loc(i,2)), 1:length(Locs));
            
            
            myCluster.Length=myCluster.AssemblyEnd-myCluster.AssemblyStart;
            
            myCluster.Insert=[C.Date{mxidx(myCluster.Loc(:,1))}]>[C.Date{mnidx(myCluster.Loc(:,1))}];
           
            myCluster.MaxCov=arrayfun(@(i) mean(mx(Locs{i})), 1:length(Locs));
            myCluster.MinCov=arrayfun(@(i) mean(mn(Locs{i})), 1:length(Locs));
            
            %isplasmid
            IsPlasmid=arrayfun(@(i) any(contains(lower(myCluster.Description{i}), 'plasmid') | startsWith(myCluster.Genes{i}, plasmid_genes)) , 1:length(myCluster.Genes), 'UniformOutput', 0);
            IsPlasmid=[IsPlasmid{:}] | myCluster.MaxCov>3;
            
            %isphage
            IsPhage=arrayfun(@(i) any(startsWith(myCluster.Genes{i}, phage_genes)) , 1:length(myCluster.Genes), 'UniformOutput', 0);
            myCluster.IsPhage=[IsPhage{:}];
            
            %isTra
            isTra=arrayfun(@(i) any(contains(myCluster.Genes{i}, 'tra')) | ...
                any(contains(myCluster.Description{i}, 'transpos')), 1:length(myCluster.Genes), 'UniformOutput', 0);
            myCluster.IsTra=[isTra{:}];
            
            
            %isRes
            %if cluster num and gene source are the same and the gene is
            %between the start and end of the cluster, isRes=true
            linIsRes=zeros(size(myCluster.IsTra))
            myCluster.IsRes=linIsRes;
            myCluster.ResGenes=cell(size(myCluster.IsRes));
            myCluster.ResPh=cell(size(myCluster.IsRes));
            
            if isfield(myCase, 'Res')
                IsRes= myCluster.Assembly==myCase.Res.Num & ...
                    myCluster.AssemblyStart<myCase.Res.Start & ...
                    myCluster.AssemblyEnd>myCase.Res.End;
                linIsRes=any(IsRes,1);
                myCluster.IsRes=linIsRes;
               myCluster.ResGenes(linIsRes)=arrayfun(@(i) ...
                myCase.Res.Gene(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
               myCluster.ResPh(linIsRes)=arrayfun(@(i) ...
                 myCase.Res.PredictedPhenotype(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
            end
            
            linIsPlas=zeros(size(myCluster.IsTra))
            myCluster.IsPlasmid=linIsRes;
            myCluster.Plasmids=cell(size(myCluster.IsPlasmid));

            if isfield(myCase, 'Plas')
                IsPlas= myCluster.Assembly==myCase.Plas.Num & ...
                    myCluster.AssemblyStart<myCase.Plas.Start & ...
                    myCluster.AssemblyEnd>myCase.Plas.End;
                linIsPlas=any(IsPlas,1);
                myCluster.IsPlasmid=linIsPlas;
               myCluster.Plasmids(linIsPlas)=arrayfun(@(i) ...
                myCase.Plas.Plasmid(IsPlas(:,i)), find(linIsPlas), 'uniformoutput',0); 
            end
            myCluster.IsPlasmid=myCluster.IsPlasmid | IsPlasmid;
            
            
            
            %number of actual definitely independent clusters
            %if the assembly is the same as the assembly of th egenes at
            %both sides
            sL=myCluster.Assembly == myCase.AssemblyNum(myCluster.Loc(:,1)-1)';
            sR=myCluster.Assembly == myCase.AssemblyNum(min(myCluster.Loc(:,2)+1, length(myCase.AssemblyNum)))';
            myCluster.N= sum(sR & sL);
            if myCluster.N<numel(myCluster.Assembly)
                myCluster.N=myCluster.N+1;
            end
            
            
            save([Path.Clusters dl 'tree' num2str(C.Num) '_' id], 'myCluster' )
            ClusterCase=[ClusterCase myCluster];
        end
        
    end
    save([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase');
end






