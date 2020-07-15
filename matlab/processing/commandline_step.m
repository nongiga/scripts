for i=1:h
    %run roary
    roary_command=[pipevar.bash_folder{i} 'roary.sh ' Path.Alignment ' ' num2str(pipevar.mm(i))];
    system(roary_command, '-echo');

    %split sequence
    split_sequence_command=sprintf("%s/split_sequence.sh %s %s %d %d", ...
        pipevar.bash_folder{i}, Path.Alignment, pipevar.read_folder{i}, pipevar.bp(i),  pipevar.mm(i));
    system(split_sequence_command, '-echo');

    %build bowtie database
    bowtie_build_command=[pipevar.bash_folder{i} 'bowtie_build.sh ' Path.Alignment ' ' num2str(pipevar.mm(i))];
    system(bowtie_build_command, '-echo');


    %run bowtie for short alignment
    if pipevar.report_multi(i)
            bowtie_align_command=[pipevar.bash_folder{i} 'bowtie_multialign.sh ' ...
            Path.Alignment ' ' num2str(pipevar.bp(i)) ' ' num2str(pipevar.mm(i))];
    else
        %read -r maindir bp threads <<<$(echo "$1 $2 $3")
        bowtie_align_command=[pipevar.bash_folder{i} 'bowtie_align.sh ' ...
            Path.Alignment ' ' num2str(pipevar.bp(i)) ' ' num2str(pipevar.mm(i))];
    end
    system(bowtie_align_command,'-echo');
end