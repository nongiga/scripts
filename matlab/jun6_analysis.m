function jun6_analysis
close all
clear all
dl=filesep;
load('gap_data', 'moptions', 'pipevar', 'Path','IsolatesNames');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

ins=1;
m=moptions{pipevar.report_multi(ins)+1};
id=[ num2str(pipevar.bp(ins)) m  ];  
load([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase')
load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')


%are there more insertions or deletions?
Inserts=arrayfun(@(c) mean(c.Insert), ClusterCase);
sum(Inserts==1)
sum(Inserts==0)
sum(Inserts>0 & Inserts<1)
%interestingly more deletions

fig=figure(1);
for iCase=1:numel(ClusterCase)
    iCase
    ax1=axes(fig)
    
    Cl=ClusterCase(iCase);
    C=Case(Cl.CaseNum);

    load([Path.Reports dl 'tree' num2str(C.Num) '_20.mat'],'myCase' )

    [mn, mnidx]=min(C.NCov, [],2);
    [mx, mxidx]=max(C.NCov, [], 2);
    
    

    
    del=arrayfun(@(i) Cl.Loc(i,1):Cl.Loc(i,2), 1:size(Cl.Loc,1), 'UniformOutput',false);
    del=[del{:}];
    
    genes=vertcat(Cl.Genes{:});
    ass_num=cellstr(num2str(myCase.AssemblyNum(del)));
    genes=join([genes, ass_num]);
    subplot(2,1,1)
    hold on
    mn=mn';mx=mx';
    mn(2,:)=nan;
    mx(2,:)=nan;
    
    H = plot(mn, mx, 'g.');
    
    arrayfun(@(i) set(H(i), 'Color', [0 0 1]), del);
    
    %H2 = plot(mn(del), mx(del), 'b.');
    title(sprintf('Num of del = %d/ %d', length(del), numel(mn)))
    
    set(H,'ButtonDownFcn',@upon_click)
    

%     subplot(2,1,2)
%     hold on
%     genesource=myCase.GeneSource;
%     
%     pos=1:numel(mn);
%     H3= plot(pos, mn./mx, 'g.');
%     H4=plot(del, mn(del)./mx(del),'b.');
%     
%     set(H3,'ButtonDownFcn',@upon_click)
%     set(H4,'ButtonDownFcn',@upon_click)
% 
% 
%     set(gca,'Xtick',pos,'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
%     set(gca,'TickLength',[0.001, 0.01])
    pause;
    clf
end


%     
function upon_click(s,~)
    %disp(s.UserData) %
    i=s.SeriesIndex;
    disp([mn(1,i) mx(1,i)])
    %plot_alignment_report(allPhages(m), hMin, hMax, t,ax2);

end


end