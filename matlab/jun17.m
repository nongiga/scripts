
close all
clear all
dl=filesep;
load('path')
AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

CD=[Path.Main dl 'deletion_reports'];


ClusterCase=[];
ErronousCase=[];

for iCase = 1: numel(Case)
    
	C=Case(iCase);
    
    load([AD dl 'tree' num2str(C.Num) '_20bp'] )
	[mn, mnidx]=min(C.Cov, [],2);
	[mx, mxidx]=max(C.Cov, [], 2);
	del=mx>0.5 & mn<0.05;
	sum(del)
	if sum(del) && length(del)<=7000 
%         figure(iCase);
% 		subplot(2,1,1)
% 		hold on
% 		plot(mn, mx, 'g.');
% 		plot(mn(del), mx(del), 'r.')
% 		title(sprintf('Num of del = %d/ %d', sum(del), numel(del)))
% 		subplot(2,1,2)
% 		hold on
% 		
% 		genesource=myCase.GeneSource;
% 		pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
% 		plot(pos, mn./mx, 'g.')
% 		plot(pos(del), mn(del)./mx(del),'r.')
%         set(gca,'Xtick',1:numel(del),'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
        
        
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
        
        save([CD dl 'tree' num2str(C.Num) '_20bp'], 'myCluster' )
        %time deletions: before or after?
        
        ClusterCase=[ClusterCase myCluster];
    elseif length(del)>7000
        disp(C.Num);
        ErronousCase=[ErronousCase C];
	end
end

save([CD  dl 'erronous_cases.mat'], 'ErronousCase');
save([CD  dl 'all_alignments20bp_cluster_report.mat'], 'ClusterCase');




