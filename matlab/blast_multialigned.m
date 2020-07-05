clear all
close all


dl=filesep;


load('../path')
load([Path.Main dl 'alignmentReports' dl 'all_alignments20bp_report.mat'])

iCase=1;
%for iCase=1:numel(Case)
C=Case(iCase);
%get all the no cov genes
zcov=all(C.NCov==0,2);

CaseFolder=[Path.Alignment dl 'tree' num2str(C.Num) dl] ;
load([Path.Main dl 'alignmentReports' dl  'tree' num2str(C.Num) '_20bp']);
myCase.GeneName(zcov);

%get these genes from the pangenome
FR=fastaread([ CaseFolder dl 'roary_output' dl 'pan_genome_reference.fa']);
FR=FR(endsWith({FR.Header}, myCase.GeneName(zcov)));

unfastq=[Path.Main dl 'blast_multialignments' dl 'tree' num2str(C.Num) '_unaligned_genes.fastq'];
fastawrite(unfastq, FR);
blastformat('Inputdb', unfastq, 'protein', false);

%get the multialigned.fastq file
FN=join([repmat({'Sample_Maccabi_Ecoli_SeqPlate'},numel(C.SeqPlate),1), C.SeqPlate'],'');

MAFile=join([repmat({CaseFolder},numel(FN),1) FN repmat({'/multialigned20.fastq'}, numel(FN), 1)],'');
system(['seqkit fq2fa ' MAFile{1} ' >tree' C.SeqPlate{1} '.fa']);

warning('off','bioinfo:blastreadlocal:NoHits');

results=blastlocal('InputQuery',  ['tree' C.SeqPlate{1} '.fa'], 'Program', 'blastn', 'Database', unfastq);
results=results(arrayfun(@(i) ~isempty(results(i).Hits), 1:numel(results)));
a=[];
for i=1:numel(results)
    a=[a, {results(i).Hits.Name}];
end
a=unique(a)

save('results', 'results', 'a')
%end