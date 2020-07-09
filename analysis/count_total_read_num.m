
dl=filesep;

roary_path = '/home/kishonystud/kishonyserver/noga/MaccabiUTI/roary' ;

tree_folders=dir([roary_path dl 'tree*']);

load('path')
AD=[ 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')


C=Case(i);
j=1;

Case.TotalReads=zeros(1, numel
[~,nreads]=system(['zcat ' roary_path dl 'tree' num2str(C.Num) dl  'Sample_Maccabi_Ecoli_SeqPlate' ...
    C.SeqPlate{j} dl 'R1_spl20.combined.trimmed.fastq.gz | wc -l'])
Case.TotalReads(j)=str2num(nreads)