
clear all; close all;
dl=filesep;

roary_path = '/home/kishonystud/kishonyserver/noga/MaccabiUTI/roary' ;

tree_folders=dir([roary_path dl 'tree*']);



for t=tree_folders'
    
    % get gene names from roary_output/pan_genome_reference.fa
    Fasta_Info=fastaread([roary_path dl t.name dl 'roary_output/pan_genome_reference.fa'])
    Gene_Info=[];
    %get info into Gene_Info
    sp=split({Fasta_Info.Header}', " ");
    Gene_Info.Sequence={Fasta_Info.Sequence}';
    Gene_Info.Source=sp(:,1);
    Gene_Info.Gene=sp(:,2);
    Gene_Info=struct2table(Gene_Info);
    
    % get all hits.ls files in directory
    Samples=dir([roary_path dl t.name dl 'Sample_Maccabi*' ]);
    Samples=Samples([Samples.isdir]);
    
    M=Gene_Info;
    
    %read hits.ls into table
    for name={Samples.name}
        name=char(name)
        T=readtable([roary_path dl t.name dl Samples(1).name dl 'hits.ls' ], ...
            'Filetype','text', 'Delimiter', 'tab', 'ReadVariableNames', false);
        HitsT=cell2table(split(table2cell(T), ' '), 'VariableNames', {name, 'Source'});
        HitsT.(name)=str2double(HitsT.(name));
        HitsT=HitsT(contains(HitsT.Source, 'SeqPlate'),:);

        M=outerjoin(M, HitsT, 'Keys', 'Source', 'MergeKeys', 1);
    end
    
end