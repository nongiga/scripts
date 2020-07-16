%loop through each pangenome
dl=filesep;
for ins=1:h
    
    SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i=1:height(SS)
        %get entries of each gene name in pangenome
        if ~exist([Path.Alignment dl TreeName dl 'pan_genome_reference.faa'], 'file')
        TreeName=['tree' SS.Case{i}]
        FastaName=[Path.Alignment dl TreeName dl 'roary_output' dl 'pan_genome_reference.fa'];
        [pangenes,~]=fastaread(FastaName);
        pangenes_s=extractBefore(pangenes, ' ');
        %loop through fasta files
        %get all fasta files with this entry name
        Header=[]; Sequence=[];
        sp=SS.seqPlates{i};
        for j=1:numel(sp)
            FileName=[Path.Alignment dl TreeName dl GlobalName sp{j} '.faa'];
            [Htemp, Stemp]=fastaread(FileName);
            Header=[Header extractBefore(Htemp, ' ')];
            Sequence=[Sequence Stemp];
        end
        
        
        [genenames_s, ia, ib]=intersect(Header,pangenes_s);
        fastawrite([Path.Alignment dl TreeName dl 'pan_genome_reference.faa'], ...
            pangenes(ib), Sequence(ia));
        end
        

    end
end