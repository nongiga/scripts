dl=filesep;
close all
load('path')
AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

CD=[Path.Main dl 'deletion_reports'];
mkdir(CD)

for iCase = 1: numel(Case)
	C=Case(iCase);
	[mn, mnidx]=min(C.Cov, [],2);
	[mx, mxidx]=max(C.Cov, [], 2);
	del=mx>0.5 & mn<0.05;
	sum(del)
	if sum(del) && ~exist([CD dl 'tree' num2str(C.Num) '_20bp'], 'file')
        
        load([AD dl 'tree' num2str(C.Num) '20bp'] )
        myCase.GeneIndel=del;
        myCase.GeneInsertion=[myCase.Date{mxidx(del)}]>[myCase.Date{mnidx(del)}];
        
        %define the beginning & end of clusters (inclusive)
        cluster_loc=find(arrayfun(@(i) del(i)~=del(i-1), 2:numel(del)));
        if mod(size(cluster_loc, 2),2)~=0
            cluster_loc=[cluster_loc length(del)];
        end
        myCase.Clusters=reshape(cluster_loc, 2, [])';
        myCase.Clusters(:,1)= myCase.Clusters(:,1)+1;
        
        %define if it is a plasmid
        plas=mx>10 & mn<0.05;
        if sum(plas)
            myCase.PlasIndel=plas;
            plas_loc=find(arrayfun(@(i) plas(i)~=plas(i-1), 2:numel(plas)))
            if mod(size(cluster_loc, 2),2)~=0
                plas_loc=[plas_loc length(del)];
            end
            
            myCase.PlasCluster=reshape(plas_loc, 2, [])';
        end
        save([CD dl 'tree' num2str(C.Num) '_20bp'] )
        %time deletions: before or after?
        
        
% 		figure(iCase);
% 		subplot(2,1,1)
% 		hold on
% 		plot(mn, mx, 'g.');
% 		plot(mn(del), mx(del), 'r.')
% 		title(sprintf('Num of del = %d/ %d', sum(del), numel(del)))
% 		subplot(2,1,2)
% 		hold on
		
% 		genesource=myCase.GeneSource;
% 		pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
% 		plot(pos, mn./mx, 'g.')
% 		plot(pos(del), mn(del)./mx(del),'r.')
%         set(gca,'Xtick',1:numel(del),'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
	end
end







