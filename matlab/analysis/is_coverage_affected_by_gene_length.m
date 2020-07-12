dl=filesep;
close all
load('path')

AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

CD=[Path.Main dl 'deletion_reports'];
mkdir(CD)
iCase=1;
%for iCase = 1: numel(Case)
C=Case(iCase);

load([AD dl 'tree' num2str(C.Num) '20bp'] )
myCase.GeneIndel=del;

%end