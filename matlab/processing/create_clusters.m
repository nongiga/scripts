
% close all
% clear all
dl=filesep;
load('gap_data');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

if ~exist('Case', 'var')
    load([Path.Reports dl 'all_alignments20_report.mat'])
end

disp("create clusters")
h=1;
phage_genes={'bet','Alpha','S_','N_','C_','4_','3_','gin','cox',...
    'B_','Q_','O_','L_','X_','Y_','R_','V_','repC','bor','ant','Alpha',...
    'ninB','cre','pol', 'cI','C2', 'MOD',...
    'P_','cII','gam','bet','exo', 'yokD', 'gag'};

plasmid_genes={'mbe', 'rop','traA', 'traB', 'traE', 'traC', 'traF', ...
    'traG', 'traH', 'traK', 'traL', 'traQ', 'traU', 'traV', 'traW'};

Path.Alignment='/home/kishonystud/kishonyserver/noga/MaccabiUTI/gap_run_prokka_reannontation/Cases';
for ins=1:h
    ClusterCase=[];
    m=moptions{pipevar.report_multi(ins)+1};
    id='20';%[ num2str(pipevar.bp(ins)) m  ];   

%     load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')        
            n_finds=0;
            n_potentials=0;
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
            
            myCluster.NumGenes=myCluster.Loc(:,2)-myCluster.Loc(:,1)+1;
            Locs=arrayfun(@(i) myCluster.Loc(i,1):myCluster.Loc(i,2), 1:length(myCluster.NumGenes), 'UniformOutput', 0);
            
            myCluster.Genes=arrayfun(@(i) C.GeneName(Locs{i}), 1:length(Locs), 'UniformOutput',false);
            myCluster.Description=arrayfun(@(i) myCase.GeneDescription(Locs{i}), 1:length(Locs), 'UniformOutput',false)';
            
            myCluster.Source=arrayfun(@(i) myCase.GeneSource(myCluster.Loc(i,1)), 1:length(Locs));
            
            myCluster.Scaffold=arrayfun(@(i) unique(myCase.AssemblyNum(myCluster.Loc(i,1):myCluster.Loc(i,2))), 1:length(Locs));
            
            myCluster.ScaffoldStart=arrayfun(@(i) myCase.GeneStart(myCluster.Loc(i,1)), 1:length(Locs));
            myCluster.ScaffoldEnd=arrayfun(@(i) myCase.GeneEnd(myCluster.Loc(i,2)), 1:length(Locs));
            
            
            myCluster.Length=myCluster.ScaffoldEnd-myCluster.ScaffoldStart;
            
            myCluster.MinDate=[C.Date{mnidx(myCluster.Loc(:,1))}];
            myCluster.MaxDate=[C.Date{mxidx(myCluster.Loc(:,1))}];
            myCluster.Insert=myCluster.MinDate<myCluster.MaxDate;
           
            myCluster.MaxCov=arrayfun(@(i) mean(mx(Locs{i})), 1:length(Locs));
            myCluster.MinCov=arrayfun(@(i) mean(mn(Locs{i})), 1:length(Locs));
            
            myCluster.Seqplate=myCluster.Source;
            %get prokka detail
            %% get info from .gff file to get location of strand on assembly

            if numel(C.SeqPlate)>2
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
                end
            end
            

            %isplasmid
            IsPlasmid=arrayfun(@(i) any(contains(lower(myCluster.Description{i}), 'plasmid') | startsWith(myCluster.Genes{i}, plasmid_genes)) , 1:length(myCluster.Genes), 'UniformOutput', 0);
            IsPlasmid=[IsPlasmid{:}] | myCluster.MaxCov>3;
            
            %isphage
            IsPhage=arrayfun(@(i) any(startsWith(myCluster.Genes{i}, phage_genes)) , 1:length(myCluster.Genes), 'UniformOutput', 0);
            myCluster.IsPhage=[IsPhage{:}];
            
            %isTra
            isTra=arrayfun(@(i) any(contains(myCluster.Description{i}, {'transpos', 'IS'})),...
                1:length(myCluster.Genes), 'UniformOutput', 0);
            myCluster.IsTra=[isTra{:}];
            
            
            %isRes
            %if cluster num and gene source are the same and the gene is
            %between the start and end of the cluster, isRes=true
            linIsRes=zeros(size(myCluster.IsTra));
            myCluster.IsRes=linIsRes;
            myCluster.ResGenes=cell(size(myCluster.IsRes));
            myCluster.ResPh=cell(size(myCluster.IsRes));
            
            if isfield(myCase, 'Res')
                IsRes= myCluster.Scaffold==myCase.Res.Num & ...
                    myCluster.ScaffoldStart<=myCase.Res.Start+50 & ...
                    myCluster.ScaffoldEnd>=myCase.Res.End-50;
                linIsRes=any(IsRes,1);
                myCluster.IsRes=linIsRes;
               myCluster.ResGenes(linIsRes)=arrayfun(@(i) ...
                myCase.Res.Gene(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
               myCluster.ResPh(linIsRes)=arrayfun(@(i) ...
                 myCase.Res.PredictedPhenotype(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
            end

            
            linIsPlas=zeros(size(myCluster.IsTra));
            myCluster.IsPlasmid=linIsRes;
            myCluster.Plasmids=cell(size(myCluster.IsPlasmid));
            
             if isfield(myCase, 'Plas') && ~isfield(myCase.Plas, 'Num')
                myCase=rmfield(myCase, 'Plas');
                save([Path.Reports dl 'tree' num2str(C.Num) '_' id  '.mat'], 'myCase')

            end           

            if isfield(myCase, 'Plas')
                sorted=sort([myCase.Plas.Start myCase.Plas.End],2);
                myCase.Plas.End=sorted(:,2); myCase.Plas.Start=sorted(:,1);
                isIntersect=~( myCluster.ScaffoldStart>myCase.Plas.End | myCluster.ScaffoldEnd<myCase.Plas.Start);
                toCalc=isIntersect==-1;
                
                
                os=max(myCluster.ScaffoldStart,myCase.Plas.Start);
                oe=min(myCluster.ScaffoldEnd,myCase.Plas.End);
                overlap=oe-os;
                overlap(~toCalc)=0;
                
                overlap=overlap./min(myCluster.Length, myCase.Plas.End-myCase.Plas.Start);
                %if there is a plasmid gene on the scaffold then the scaffold is
                %a plasmid

                IsPlas= myCluster.Scaffold==myCase.Plas.Num;
                linIsPlas=any(IsPlas,1);
                myCluster.IsPlasmid=linIsPlas;
                myCluster.Plasmids(linIsPlas)=arrayfun(@(i) ...
                    myCase.Plas.Plasmid(IsPlas(:,i)), find(linIsPlas), 'uniformoutput',0); 
            end
            myCluster.IsPlasmid=myCluster.IsPlasmid | IsPlasmid;
            
            
            
            %number of actual definitely independent clusters
            %if the assembly is the same as the assembly of th egenes at
            %both sides
            sL=myCluster.Scaffold == myCase.AssemblyNum(myCluster.Loc(:,1)-1)';
            sR=myCluster.Scaffold == myCase.AssemblyNum(min(myCluster.Loc(:,2)+1, length(myCase.AssemblyNum)))';
            myCluster.Fragmented=~(sR & sL);
            myCluster.N= sum(sR & sL);
            if myCluster.N<numel(myCluster.Scaffold)
                myCluster.N=myCluster.N+1;
            end
            
            
            save([Path.Clusters dl 'tree' num2str(C.Num) '_' id], 'myCluster' )
            ClusterCase=[ClusterCase myCluster];
        end
        
    end
    save([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase');
end






