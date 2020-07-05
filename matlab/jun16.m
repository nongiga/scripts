dl=filesep;
close all
load('path')
AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

for iCase = 1: 10
	C=Case(iCase);
	mn=min(C.Cov, [],2);
	mx=max(C.Cov, [], 2);
	del=mx>0.5 & mn<0.05;
	sum(del)
	if sum(del)
		figure(iCase);
		subplot(2,1,1)
		hold on
		plot(mn, mx, 'g.');
		plot(mn(del), mx(del), 'r.')
		title(sprintf('Num of del = %d/ %d', sum(del), numel(del)))
		subplot(2,1,2)
		hold on
		load([AD dl 'tree' num2str(C.Num) '20bp'] )
		genesource=myCase.GeneSource;
		pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
		plot(pos, mn./mx, 'g.')
		plot(pos(del), mn(del)./mx(del),'r.')
        set(gca,'Xtick',1:numel(del),'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
	end
end







