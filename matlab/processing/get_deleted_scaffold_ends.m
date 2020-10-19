load('cluster_table')
load('gap_data.mat')
dl=filesep;
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';
Path.Alignment='Cases';
id='20';
load([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase')
if ~exist('Case', 'var')
    load([Path.Reports dl 'all_alignments20_report.mat'])
end

fullscaf=struct('Header',{},'Sequence',{});
bp20scaf=struct('Header',{},'Sequence',{});
k=1;
for i=1:numel(ClusterCase)
    myCluster=ClusterCase(i);
    C=Case( ClusterCase(i).CaseNum);
    
     FolderName=myCluster.Source{1}(1:end-6);
    IsolateDir=['assembly' dl FolderName dl 'contigs.fasta' ];
    try
        FASTAData=fastaread(IsolateDir);
    catch
        continue
    end
    for j=1:numel(myCluster.Scaffold)
        FolderName=myCluster.Source{j}(1:end-6);
        if ~contains(IsolateDir, FolderName)
            IsolateDir=['assembly' dl FolderName dl 'contigs.fasta' ];
            FASTAData=fastaread(IsolateDir);
        end

        fastaloc=find(startsWith({FASTAData.Header}, ['NODE_' num2str(myCluster.Scaffold(j)) '_']));
        
        header=[FolderName '-' FASTAData(fastaloc).Header];
        
        fullscaf(k).Header=header; 
        fullscaf(k).Sequence=FASTAData(fastaloc).Sequence;
        
        bp20scaf(2*k-1).Header=[header '-START'];
        bp20scaf(2*k-1).Sequence=FASTAData(fastaloc).Sequence(1:20);
        
        bp20scaf(2*k).Header=[header '-END'];
        bp20scaf(2*k).Sequence=FASTAData(fastaloc).Sequence(end-20:end);
        
        k=k+1;
    end
end
fastawrite('cluster_scaffold_ends.fasta', bp20scaf);
fastawrite('cluster_scaffolds.fasta', fullscaf);



% zcat ../Mathews_processing/Filtered_data/Sample_Maccabi_Ecoli_SeqPlate10_D5/R1_combined.trimmed.fastq.gz | seqkit grep -s -p 'GGGACCGCGGTCCCACTCGT'
