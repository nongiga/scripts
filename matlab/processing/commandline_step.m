for ins=1:h
    %run prokka
    AlPath=[pwd dl Path.Alignment];
    prokka_command=[pipevar.bash_folder{ins} 'prokka.sh ' pwd dl Path.Assembly ' ' Path.Prokka ' ' num2str(pipevar.mm(ins))];
    system(prokka_command, '-echo');
    
    initialize_case_folders
    
    %run roary
    roary_command=[pipevar.bash_folder{ins} 'roary.sh ' AlPath ' ' num2str(pipevar.mm(ins))];
    disp(roary_command)
    system(roary_command, '-echo');

    %split sequence
    split_sequence_command=sprintf("%s/split_sequence.sh %s %s %d %d", ...
        pipevar.bash_folder{ins}, AlPath, pipevar.read_folder{ins}, pipevar.bp(ins),  pipevar.mm(ins));
    disp(split_sequence_command);
    system(split_sequence_command, '-echo');

    %build bowtie database
    bowtie_build_command=[pipevar.bash_folder{ins} 'bowtie_build.sh ' AlPath ' ' num2str(pipevar.mm(ins))];
    system(bowtie_build_command, '-echo');


    %run bowtie for short alignment
    if pipevar.report_multi(ins)
            bowtie_align_command=[pipevar.bash_folder{ins} 'bowtie_multialign.sh ' ...
            AlPath ' ' num2str(pipevar.bp(ins)) ' ' num2str(pipevar.mm(ins))];
    else
        %read -r maindir bp threads <<<$(echo "$1 $2 $3")
        bowtie_align_command=[pipevar.bash_folder{ins} 'bowtie_align.sh ' ...
            AlPath ' ' num2str(pipevar.bp(ins)) ' ' num2str(pipevar.mm(ins))];
    end
     system(bowtie_align_command,'-echo');
end