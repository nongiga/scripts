close all
clear all
dl=filesep;
load('gap_data');
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
figure(1)
for iCase=1:numel(ClusterCase)
    iCase
    
    Cl=ClusterCase(iCase);
    C=Case(Cl.CaseNum);

    load([Path.Reports dl 'tree' num2str(C.Num) '_20.mat'] )

    [mn, mnidx]=min(C.NCov, [],2);
    [mx, mxidx]=max(C.NCov, [], 2);
    
    del=arrayfun(@(i) Cl.Loc(i,1):Cl.Loc(i,2), 1:size(Cl.Loc,1), 'UniformOutput',false);
    del=[del{:}];
    
    genes=vertcat(Cl.Genes{:});
    ass_num=cellstr(num2str(myCase.AssemblyNum(del)));
    genes=join([genes, ass_num]);
    subplot(2,1,1)
    hold on
    plot(mn, mx, 'g.');
    plot(mn(del), mx(del), 'b.')
    title(sprintf('Num of del = %d/ %d', length(del), numel(mn)))

    subplot(2,1,2)
    hold on
    genesource=myCase.GeneSource;
%     lastdot_pos = cellfun(@(gs) find(gs == '_', 1, 'last'), genesource);
%     unique(arrayfun(@(i) genesource{i}(1:lastdot_pos(i)-1), 1:length(lastdot_pos), 'UniformOutput', false));
%    pos=cellfun(@(x) str2num(extractAfter(x(end-6:end), '_')),genesource);
    pos=1:numel(mn);
    plot(pos, mn./mx, 'g.')
    plot(del, mn(del)./mx(del),'b.')
    

    set(gca,'Xtick',pos,'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
    set(gca,'TickLength',[0.001, 0.01])
    pause;
    clf
end

