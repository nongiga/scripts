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

%% get all gene names

all_genes=[ClusterCase.Genes];
all_descriptions=vertcat(ClusterCase.Description);
all_genes=vertcat(all_genes{:});
all_descriptions=vertcat(all_descriptions{:});

unidentified_genes=all_genes(contains(all_descriptions, 'hypothetical protein'));

all_sequences=[];

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];

casenums=arrayfun(@(i) repelem({ClusterCase(i).Num},cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];
all_locs=[];
%% get sequence of all genes
if ~exist('clustered_genes_table.mat', 'file')
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
else
    load('clustered_genes_table.mat', 'all_sequences')
end


%% make fasta files and run cd-hit
entry_name=join( [all_genes, casenums'], '-') ;
mkdir('blast');
Path.Blast='blast';
cd(Path.Blast)
blastfa='undefined.fa';
p_seq=nt2aa(all_sequences);
blastfaa='undefined_p.fa';

if ~exist('undefined_p.fa.groups', 'file')
    fastawrite(blastfaa, entry_name, p_seq);
    system('iterative_cdhit -m undefined_p.fa -f undefined_clustered_filtered.fa -c undefined_clustered -v');
end

%% process cd-hit results to cluster genes

ClusterT=readtable('undefined_p.fa.groups', 'FileType', 'text');
hypothetical_loc=ismember(all_descriptions, 'hypothetical protein');
clustered_descriptions=all_descriptions;
for i=1:height(ClusterT)
    gene_loc=ismember(entry_name, ClusterT{i,:});
    
    c= ismember(all_descriptions(gene_loc), 'hypothetical protein');
    % if all the proteins are hypothetical define them as cluster
    if all(c)
        clustered_descriptions(gene_loc)={['cluster ' num2str(i)]};
    %if one or more of the genes is defined but the rest are not
    elseif sum(c)
        % apply the description to all hypothetical proteins
        %divide descriptions to 
        descriptions=all_descriptions(gene_loc & ~hypothetical_loc);
        disp(descriptions)
        clustered_descriptions(gene_loc & hypothetical_loc)=descriptions(1);
        
    end
    
end

%% now all the rest of the undefined genes that were not clustered will be named based
% on their source

[u, uidx]=unique(clustered_descriptions);

table_gene_names=all_genes;
table_gene_names(contains(all_genes, 'group'))=entry_name(contains(all_genes, 'group'));


% find those genes with the same name but different description and add _d
% to them to differentiate
unique_gene_names=table_gene_names(uidx);
[~, ind] = unique(unique_gene_names);
duplicate_ind = setdiff(1:size(unique_gene_names, 1), ind);
table_gene_names(uidx(duplicate_ind))=strcat(unique_gene_names(duplicate_ind), '_d');

for i=1:numel(u)
    %find location of description
    %name all table_gene_names according to the same table_gene_name with
    %this description
    
    
end


T=table('Size', [numel({ClusterCase.Num}) numel(u)], 'VariableNames', table_gene_names(uidx), ...
    'VariableType', repelem({'uint8'}, numel(u)));
T.Properties.VariableDescriptions=u;
T.Properties.RowNames={ClusterCase.Num};


is_plasmid=arrayfun(@(i) arrayfun(@(j) repmat(ClusterCase(i).IsPlasmid(j), ClusterCase(i).GeneNum(j),1), 1:numel(ClusterCase(i).GeneNum), 'UniformOutput', false), 1:numel(ClusterCase), 'UniformOutput', false);
is_plasmid=[is_plasmid{:}];
is_plasmid=vertcat(is_plasmid{:});


for i=1:numel(casenums)
    T{casenums{i}, table_gene_names{i}}=uint8(is_plasmid(i))+1;
end

save('clustered_genes_table.mat','T', 'all_sequences');

%% simulate random permutation of numbers in vector
seed=100;
Tmat=table2array(T);
for i=1:height(T)
    seed=seed+1;
    s=rng(seed);
    Tmat(i,:)=Tmat(i, randperm(width(T)));
end


binRange=0:1:10;

sim=histcounts(sum(Tmat>0), [binRange Inf]);
real=histcounts(sum(T{:,:}>0),[binRange Inf]);

bar(binRange, [sim;real]')

all_descriptions(uidx(sum(T{:,:}>0)>8))
