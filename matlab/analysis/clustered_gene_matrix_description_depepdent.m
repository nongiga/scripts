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
load([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')

%% get all gene names

all_genes=[ClusterCase.Genes];
all_descriptions=vertcat(ClusterCase.Description);
all_genes=vertcat(all_genes{:});
all_descriptions=vertcat(all_descriptions{:});



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


%% cd-hit
entry_name=join( [all_genes, casenums'], '-') ;
mkdir('blast');
Path.Blast='blast';
cd(Path.Blast)
blastfa='undefined.fa';
p_seq=nt2aa(all_sequences);
blastfaa='undefined_p.fa';
fastawrite(blastfaa, entry_name, p_seq);

system('iterative_cdhit -m undefined_p.fa -f undefined_clustered_filtered.fa -c undefined_clustered -v');


%% process cd-hi results to cluster genes

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
        d_g=d_g{1};
        last_pos = find(d_g == '-', 1, 'last')
        name=d_g(1:last_pos-1);
        clustered_genes(ismember(entry_name, genes(c)))={name};
        
    end
    
end

% now all the rest of the undefined genes that were not clustered will be named based
% on their source
still_undefined=contains(clustered_genes, 'group')
clustered_genes(still_undefined)=entry_name(still_undefined);


clean_genes=extractBefore(strcat(all_genes, '_'), '_');
short_names=cellfun(@length, clean_genes)<2;

% if clean_gene is the same name and dscription is the same
%then give the same name
clean_genes(short_names)=all_genes(short_names);


[u, uidx]=unique(clean_genes);


T=table('Size', [numel({ClusterCase.Num}) numel(u)], 'VariableNames', u, ...
    'VariableType', repelem({'uint8'}, numel(u)));
T.Properties.VariableDescriptions=all_descriptions(uidx);
T.Properties.RowNames={ClusterCase.Num};

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];


casenums=arrayfun(@(i) repelem({ClusterCase(i).Num},cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];

is_plasmid=arrayfun(@(i) arrayfun(@(j) repmat(ClusterCase(i).IsPlasmid(j), ClusterCase(i).GeneNum(j),1), 1:numel(ClusterCase(i).GeneNum), 'UniformOutput', false), 1:numel(ClusterCase), 'UniformOutput', false);
is_plasmid=[is_plasmid{:}];
is_plasmid=vertcat(is_plasmid{:});


for i=1:numel(casenums)
    T{casenums{i}, clean_genes{i}}=uint8(is_plasmid(i))+1;
end

save('clustered_genes_table.mat','T');

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

all_descriptions(uidx(sum(T{:,:}>0)>7))
