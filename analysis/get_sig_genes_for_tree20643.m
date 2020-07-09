close all
clear all

dl=filesep;

AD='alignmentReports';
CD='deletion_reports';

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

C=Case(179);
C.Num
load([AD dl 'tree' num2str(C.Num) '_20bp'] )

[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);

del=mx>0.5 & mn<0.05 ;
sig=C.pval==0;


figure(1);
subplot(2,1,1)
hold on
plot(mn, mx, 'g.');
plot(mn(del), mx(del), 'b.')
plot(mn(sig), mx(sig), 'r.')

title(sprintf('Num of del = %d/ %d', sum(del), numel(del)))

subplot(2,1,2)
hold on
genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')
plot(pos(del), mn(del)./mx(del),'b.')
plot(pos(sig), mn(sig)./mx(sig),'r.')



sig_loc=find(C.pval<10^-5);
genesource(sig_loc(1))
myCase.GeneName(sig_loc(1))
myCase.Reads(sig_loc(1),:)
myCase.NCov(sig_loc(1),:)
myCase.GeneLength(sig_loc(1))

%13
 