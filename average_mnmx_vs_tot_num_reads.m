close all
clear all

dl=filesep;

AD='alignmentReports';
CD='deletion_reports';

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

iCase=1;
    
C=Case(iCase);


[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);

del=mx>0.5 & mn<0.05 ;
sig=C.pval==0;
v=arrayfun(@(c) var(c.NCov(:)), Case);

r=arrayfun(@(c) min(sum(c.Reads)), Case);
figure(1);
plot(log10(r), v, 'g.')
xlabel('log min # of reads on isolate in case');
ylabel('variance of normalized counts');


load([AD dl 'mean_vs_pval'], 'x', 'y')
 