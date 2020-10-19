
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
        del=any(C.NCov>0.5,2) & any(C.NCov<0.05,2);
        sum(del)
        if sum(del) && numel(del)<6000 && exist([Path.Reports dl 'tree' num2str(C.Num) '_' id  '.mat'], 'file')
            load([Path.Reports dl 'tree' num2str(C.Num) '_' id  '.mat'] )
            mC=myCase;
            delf=find(del);
            gs=cellfun(@(g) g(1:find(g=='_', 1,'last')-1), mC.GeneSource, 'UniformOutput', false);
            %define the beginning & end of clusters (inclusive)
            cl=get_cluster_location(C.NCov,del, mC.AssemblyNum, gs);
            
            if isempty(cl),     continue; end
           
            Locs=arrayfun(@(i) cl(i,1):cl(i,2), 1:size(cl,1), 'UniformOutput', 0);
            MinIso=arrayfun(@(i) unique(C.NCov(intersect(Locs{i},delf),:)<0.10,'rows'),1:length(Locs), 'UniformOutput',0);
            MaxIso=arrayfun(@(i) unique(C.NCov(intersect(Locs{i},delf),:)>0.10,'rows'),1:length(Locs), 'UniformOutput',0);
            
            mCl=struct('Loc', {cl}, 'CaseNum', iCase, 'Num', C.Num,...
                'NumGenes', {cl(:,2)-cl(:,1)+1}, ...
                'Genes', {arrayfun(@(i) C.GeneName(Locs{i}), 1:length(Locs), 'UniformOutput',false)}, ...
                'Description', {arrayfun(@(i) mC.GeneDescription(Locs{i}), 1:length(Locs), 'UniformOutput',false)'},...
                'Source', {arrayfun(@(i) mC.GeneSource(cl(i,1)), 1:length(Locs))},...
                'Scaffold', {arrayfun(@(i) unique(mC.AssemblyNum(Locs{i})), 1:length(Locs))},...
                'Start', {arrayfun(@(i) mC.GeneStart(Locs{i}(1)), 1:length(Locs))},...
                'End', {arrayfun(@(i) mC.GeneEnd(Locs{i}(end)), 1:length(Locs))},...
                'MinDate', {cellfun(@(mi) C.Date{mi},MinIso)},...
                'MaxDate', {cellfun(@(mi) C.Date{mi},MaxIso)},...
                'MaxCov', {arrayfun(@(i) mean(C.NCov(Locs{i},find(MaxIso{i}))), 1:length(Locs), 'UniformOutput',0)},...
                'MinCov', {arrayfun(@(i) mean(C.NCov(Locs{i},find(MinIso{i}))), 1:length(Locs), 'UniformOutput',0) }, ...
                'MinPlate', {cellfun(@(mi) C.SeqPlate{mi},MinIso, 'UniformOutput',0)}, ...
                'MaxPlate',  {cellfun(@(mi) C.SeqPlate{mi},MaxIso,'UniformOutput',0)});
            
            mCl.Length=mCl.End-mCl.Start;
            mCl.Insert=mCl.MinDate<mCl.MaxDate;

            %get prokka detail
            %% get info from .gff file to get location of strand on assembly
% 
%             if numel(C.SeqPlate)>2
%                 % check if genes found at different, longer scaffolds
%                 find_longer_scaffold(Path, myCluster, C, GlobalName)
%             end
            

            %isphage
            IsPhage=arrayfun(@(i) any(startsWith(mCl.Genes{i}, phage_genes)) , 1:length(mCl.Genes), 'UniformOutput', 0);
            mCl.IsPhage=[IsPhage{:}];
            
            %isTra
            isTra=arrayfun(@(i) any(contains(mCl.Description{i}, {'transpos', 'IS'})),...
                1:length(mCl.Genes), 'UniformOutput', 0);
            mCl.IsTra=[isTra{:}];
            
            
            %isRes
            %if cluster num and gene source are the same and the gene is
            %between the start and end of the cluster, isRes=true
%             linIsRes=zeros(size(myCluster.IsTra));
%             myCluster.IsRes=linIsRes;
%             myCluster.ResGenes=cell(size(myCluster.IsRes));
%             myCluster.ResPh=cell(size(myCluster.IsRes));
%             
%             if isfield(myCase, 'Res')
%                 IsRes= myCluster.Scaffold==myCase.Res(mnidx(cl(:,1))).Num & ...
%                     myCluster.ScaffoldStart<=myCase.Res.Start+50 & ...
%                     myCluster.ScaffoldEnd>=myCase.Res.End-50;
%                 linIsRes=any(IsRes,1);
%                 myCluster.IsRes=linIsRes;
%                myCluster.ResGenes(linIsRes)=arrayfun(@(i) ...
%                 myCase.Res.Gene(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
%                myCluster.ResPh(linIsRes)=arrayfun(@(i) ...
%                  myCase.Res.PredictedPhenotype(IsRes(:,i)), find(linIsRes), 'uniformoutput',0); 
%             end
            
            %isplasmid
            
            mCl.IsPlasmid=zeros(size(mCl.IsTra));
            mCl.Plasmids=cell(size(mCl.IsPlasmid));
            if isfield(mC, 'Plas')
                IsPlas= mCl.Scaffold==mC.Plas.Num;
                linIsPlas=any(IsPlas,1);
                mCl.IsPlasmid=linIsPlas;
                mCl.Plasmids(linIsPlas)=arrayfun(@(i) mC.Plas.Plasmid(IsPlas(:,i)), find(linIsPlas), 'uniformoutput',0); 
            end
            IsPlasmid=arrayfun(@(i) any(contains(lower(mCl.Description{i}), 'plasmid') | startsWith(mCl.Genes{i}, plasmid_genes)) , 1:length(mCl.Genes), 'UniformOutput', 0);
            mCl.IsPlasmid=mCl.IsPlasmid | [IsPlasmid{:}] | arrayfun(@(mc) max(mc{1})>3,mCl.MaxCov);

            %number of actual definitely independent clusters
            %if the assembly is the same as the assembly of th egenes at
            %both sides
            sL=mCl.Scaffold == mC.AssemblyNum(mCl.Loc(:,1)-1)';
            sR=mCl.Scaffold == mC.AssemblyNum(min(mCl.Loc(:,2)+1, length(mC.AssemblyNum)))';
            
            mCl.Fragmented=~(sR & sL);
            
            mCl.Context=zeros(size(sL));
            mCl.Context= uint8(sR) +uint8(sL);
            
            mCl.N= sum(sR & sL);
            if mCl.N<numel(mCl.Scaffold)
                mCl.N=mCl.N+1;
            end
            ClusterCase=[ClusterCase mCl];
        end
        
    end
    save([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase');
end





