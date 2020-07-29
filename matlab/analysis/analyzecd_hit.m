%% set up variables
Path.Reports='alignmentReports';
Path.Blast='blast';
dl=filesep; id='20';

cellsizes=cellfun(@sum, {ClusterCase.GeneNum}, 'UniformOutput', false);
cellsizes=[cellsizes{:}];

casenums=arrayfun(@(i) repelem(i,cellsizes(i)), 1:numel(cellsizes), 'UniformOutput', false);
casenums=[casenums{:}];
all_locs=[]; all_sequences=[];

entry_name=join( [geneTable(:,1), cellstr(num2str(casenums'))], '-') ;

%% get sequence of all genes
blastfaa=[Path.Blast dl 'deleted_sequences.fa'];
if ~exist(blastfaa, 'file')
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

    % make fasta files and run cd-hit
    
    mkdir(Path.Blast);
    p_seq=nt2aa(all_sequences);
    blastfaa=[Path.Blast dl 'deleted_sequences.fa'];
    fastawrite(blastfaa, entry_name, p_seq);
end
if ~exist('blast/deleted_sequences.fa.groups', 'file')
    
    system('iterative_cdhit -m deleted_sequences.fa -f deleted_sequences.fa -c deleted_sequences -v');
end

%% process cd-hit results to cluster genes

ClusterT=readtable('blast/undefined_p.fa.groups', 'FileType', 'text');
hypothetical_loc=ismember(geneTable(:,2), 'hypothetical protein');
clustered_descriptions=geneTable(:,2);
for i=1:height(ClusterT)
    gene_loc=ismember(entry_name, ClusterT{i,:});
    
    c= ismember(geneTable(gene_loc,2), 'hypothetical protein');
    % if all the proteins are hypothetical change desc
    if all(c)
        clustered_descriptions(gene_loc)={['cluster ' num2str(i)]};
    %if one or more of the genes is defined but the rest are not
    elseif sum(c)
        % apply the description to all hypothetical proteins
        %divide descriptions to 
        descriptions=geneTable(gene_loc & ~hypothetical_loc,2);
        disp(descriptions)
        clustered_descriptions(gene_loc & hypothetical_loc)=descriptions(1);
        
    end
    
end

%% now all the rest of the undefined genes that were not clustered will be named based
% on their source

[u, uidx]=unique(clustered_descriptions);

table_gene_names=geneTable(:,1);
table_gene_names(contains(geneTable(:,1), 'group'))=entry_name(contains(geneTable(:,1), 'group'));


% find those genes with the same name but different description and add _d
% to them to differentiate
unique_gene_names=table_gene_names(uidx);
[~, ind] = unique(unique_gene_names);
duplicate_ind = setdiff(1:size(unique_gene_names, 1), ind);
table_gene_names(uidx(duplicate_ind))=strcat(unique_gene_names(duplicate_ind), '_d');

