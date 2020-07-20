for ins=1:h
    %run prokka
    prokka_command=[pipevar.bash_folder{ins} 'prokka.sh ' pipevar.assembly_folder{ins} ' ' Path.Prokka ' ' num2str(pipevar.mm(ins))];
    system(prokka_command, '-echo');
    
    initialize_case_folders

    %run roary
    roary_command=[pipevar.bash_folder{ins} 'roary.sh ' Path.Alignment ' ' num2str(pipevar.mm(ins))];
    system(roary_command, '-echo');

    %split sequence
    split_sequence_command=sprintf("%s/split_sequence.sh %s %s %d %d", ...
        pipevar.bash_folder{ins}, Path.Alignment, pipevar.read_folder{ins}, pipevar.bp(ins),  pipevar.mm(ins));
    system(split_sequence_command, '-echo');

    %build bowtie database
    bowtie_build_command=[pipevar.bash_folder{ins} 'bowtie_build.sh ' Path.Alignment ' ' num2str(pipevar.mm(ins))];
    system(bowtie_build_command, '-echo');


    %run bowtie for short alignment
    if pipevar.report_multi(ins)
            bowtie_align_command=[pipevar.bash_folder{ins} 'bowtie_multialign.sh ' ...
            Path.Alignment ' ' num2str(pipevar.bp(ins)) ' ' num2str(pipevar.mm(ins))];
    else
        %read -r maindir bp threads <<<$(echo "$1 $2 $3")
        bowtie_align_command=[pipevar.bash_folder{ins} 'bowtie_align.sh ' ...
            Path.Alignment ' ' num2str(pipevar.bp(ins)) ' ' num2str(pipevar.mm(ins))];
    end
%     system(bowtie_align_command,'-echo');
end