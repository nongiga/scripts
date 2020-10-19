%zcat ../Mathews_processing/Filtered_data/Sample_Maccabi_Ecoli_SeqPlate10_D5/R*_combined.trimmed.fastq.gz | seqkit grep -s -p 'GGGACCGCGGTCCCACTCGT' |seqkit stats

forw=fastaread('cluster_scaffold_ends.fasta');

%find matching sequences
for i=1:numel(forw)
    sp=extractBefore(forw(i).Header, '-');
    seqf=forw(i).Sequence;
    cmdfind=['zcat ../Mathews_processing/Filtered_data/' sp '/R*_combined.trimmed.fastq.gz | seqkit locate -r -p ' seqf ' | cut -f 1 '];
    system([cmdfind '>cluster_scaffold_ends/' forw(i).Header '.txt']);
    cmdlst=['zcat ../Mathews_processing/Filtered_data/' sp '/R*_combined.trimmed.fastq.gz | seqkit grep -f cluster_scaffold_ends/' forw(i).Header '.txt '];
    system([cmdlst '>cluster_scaffold_ends/' forw(i).Header '.fastq']);
    
    
    system(['/media/kishonylab/KishonyStorage/Apps/SPAdes-3.14.1-Linux/bin/spades.py --isolate -s cluster_scaffold_ends/' forw(i).Header '.fastq -o scaffold_end_assembly/'  forw(i).Header]);
    
    % blast against the scaffolds of the assembly
    
end

%assemble matching sequences
