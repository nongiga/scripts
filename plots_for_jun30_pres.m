close all
clear all

dl=filesep;

AD='alignmentReports';
CD='deletion_reports';

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

C=Case(1);
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
plot(mn(sig), mx(sig), 'r.')
plot(mn(del), mx(del), 'b.')
title('Normalized coverage of Genes in Case #1001239')
xlabel('Min normalized coverage')
ylabel('Max normalized coverage');



subplot(2,1,2)
hold on
genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')
plot(pos(sig), mn(sig)./mx(sig),'r.')
plot(pos(del), mn(del)./mx(del),'b.')
title('Fold change from minimum to maximum normalized coverage over genome')
xlabel('Gene Location in Assembly')
ylabel('Min/Max Normalized Coverage');



sig_loc=find(C.pval==0 & del);
genesource(sig_loc(1))
myCase.GeneName(sig_loc(1))
myCase.Reads(sig_loc(1),:)
myCase.NCov(sig_loc(1),:)
myCase.GeneLength(sig_loc(1))

figure(2);
clf
plot(myCase.BaseCov{4, 1}, 'r')
hold on;
plot(myCase.BaseCov{4, 2}, 'b')
title(sprintf("Coverage on gene # %s in case # %s", myCase.GeneName{4}, num2str(myCase.Num)))

baminfo('/home/kishonystud/kishonyserver/noga/MaccabiUTI/roary/tree1001239/Sample_Maccabi_Ecoli_SeqPlate17_G10/aligned_sorted20.bam')
br=bamread('/home/kishonystud/kishonyserver/noga/MaccabiUTI/roary/tree1001239/Sample_Maccabi_Ecoli_SeqPlate17_G10/aligned_sorted20.bam',35, [1 933])
%13
 