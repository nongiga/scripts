clear all; close all;
%% set-up
dl=filesep;
load('gap_data', 'moptions', 'pipevar', 'Path');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

ins=1;
m=moptions{pipevar.report_multi(ins)+1};
id=[ num2str(pipevar.bp(ins)) m  ];  

load([Path.Clusters  dl 'all_clusters' id '_report.mat'], 'ClusterCase')
%load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')

%% get all unidentified genes

all_genes=[ClusterCase.Genes];
all_genes=vertcat(all_genes{:});
unidentified_genes=all_genes(contains(all_genes, 'group_'));

all_sequences=[];

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];

casenums=arrayfun(@(i) repelem({ClusterCase(i).Num},cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];
all_locs=[];
%% get sequence of all unidentified genes
for iCC=1:numel(ClusterCase)
    CC=ClusterCase(iCC);
    load([Path.Reports dl 'tree' CC.Num '_' id '.mat'])
    
    locs = arrayfun(@(i) CC.Loc(i,1):CC.Loc(i,2), ...
        1:size(CC.Loc, 1), 'UniformOutput', 0);
    locs=[locs{:}];
    all_locs= [all_locs locs];
    cluster_sequences=myCase.Sequence(locs);
    all_sequences=[ all_sequences; cluster_sequences];

end


%% BLAST
entry_name=join( [all_genes, casenums'], '-') ;
mkdir('blast');
Path.Blast='blast';


cd(Path.Blast)
blastfa='undefined.fa';

p_seq=nt2aa(all_sequences);

fastawrite(blastfa, entry_name, all_sequences);
% blastformat('Inputdb', blastfa, 'Protein', 0); 
% blast_results=blastlocal('InputQuery', blastfa, 'Database',blastfa, 'Program', 'blastn')


blastfaa='undefined_p.fa';
fastawrite(blastfaa, entry_name, p_seq);



%blastformat('Inputdb', blastfaa, 'Protein', 1); 
%blast_results=blastlocal('InputQuery', blastfaa, 'Database',blastfaa, 'Program', 'blastp')

%save('blast_results', 'blast_results');


%% process blast_results
% match_table={};
% for i=1:numel(blast_results)
%     br=blast_results(i);
%     matches={};
%     for j=1:numel(br.Hits)
%         if br.Hits(j).HSPs(1).Identities.Percent>98
%             matches=[matches ;{br.Hits(j).Name}];
%         end
%     end
%     match_table=[match_table; {br.Query} {matches}];
%     
% end

%%  ughh just throw it out the window and use cd-hit. why re-inventthe wheel
% just because Idan hasn't heard of it?


system('iterative_cdhit -m undefined_p.fa -f undefined_clustered_filtered.fa -c undefined_clustered -v');

ClusterT=readtable('undefined_p.fa.groups', 'FileType', 'text');


clustered_genes=all_genes;
for i=1:height(ClusterT)
    genes=ClusterT{i,:};
    genes=genes(~cellfun(@isempty, genes));
    c=contains(genes, 'group');
    % if all the proteins are hypothetical define them as cluster
    if all(c)
        clustered_genes(ismember(entry_name, genes))={['cluster_' num2str(i)]};
    %if one or more of the genes is defined but the rest are not
    elseif sum(c)
        %extract gene name from first defined group and apply it to all
        d_g=genes(~c);
        name=extractBefore(d_g(1), '-');
        clustered_genes(ismember(entry_name, genes(c)))=name;
        
    end
    
end

% now all the rest of the undefined genes that were not clustered will be named based
% on their source
still_undefined=contains(clustered_genes, 'group')
clustered_genes(still_undefined)=entry_name(still_undefined);



