for ins=1:h
    
    if ~pipevar.redo(ins),  continue; end
    
    SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i = 1:height(SS)
        %make dir
        sp=SS.seqPlates{i};
        for i=1:numel(sp)
            dirname=[Path.Assembly dl GlobalName sp{i}];
            if ~exist(dirname, 'dir'), mkdir(dirname); end
            system(['cp -u ' pipevar.assembly_folder{ins} dl '*' GlobalName sp{i} dl 'spades_assembly' ...
                dl 'contigs.fasta ' dirname]);
        end
    end
end