dl=filesep;
AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20_report.mat'], 'Case')

iCase=3;
    
C=Case(iCase);

load([AD dl 'tree' num2str(C.Num) '_20.mat'] )
[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);

del=mx>0.5 & mn<0.05 ;
sum(del)
sig=myCase.pval==0;

figure(iCase);

subplot(3,1,1)
hold on
plot(mn, mx, 'g.');
plot(mn(sig), mx(sig), 'r.')
plot(mn(del), mx(del), 'b.')
title(sprintf('Num of del = %d/ %d', sum(del), numel(del)))

subplot(3,1,2)
hold on
genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')
plot(pos(del), mn(del)./mx(del),'b.')
plot(pos(sig), mn(sig)./mx(sig),'r.')
set(gca,'Xtick',1:numel(del),'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)

genelength=myCase.GeneLength;
subplot(3,1,3)
hold on
plot(genelength, mn./mx, 'g.')
plot(genelength(del), mn(del)./mx(del),'r.')



%is it gene length that is responsible for all the false positives?
% look at gene length in dels vs no dels
figure;
subplot(2,1,1)
hold on
histogram(genelength, 0:20:1600)
title('gene lengths')

subplot(2,1,2)
histogram(genelength(del), 0:20:1600)
title('gene length of positives')

% for C.Num=1001239, there is a very dense cluster of genes w/ mn=0 mx>0.5.
% Did we put the threshold in the right place?

figure;
histogram(mx(del & mn==0), 11)

%plot fold change distribution
figure;
histogram(mn./mx)
histogram(mn(del)./mx(del), 'FaceColor', 'r')


%basemean vs foldchange
 plot(log2(mn./mx),mean(C.Cov,2), '.')
 
 