for ins=1:h
    
    if ~pipevar.redo(ins),  continue; end
    
    SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i = 1:height(SS)
        %make dir
        sp=SS.seqPlates{i};

        % if same dirs exist here and there copy over all
        %dir([pipevar.ref_db_path{1} dl 'Cases' dl 'tree' SameStrains.Case{i}])
        pastTreeDir = dir([pipevar.ref_db_path{1} dl 'Cases' dl ...
            'tree' SS.Case{i}]);
        pastTreeDir(ismember( {pastTreeDir.name}, {'.', '..'}) | [pastTreeDir.isdir]==0) = []; 

        expectedDirs=[strcat(GlobalName, sp), {'roary_output'}];

        %if intersection between lists is perfect copy over
        if all(ismember({pastTreeDir.name}, expectedDirs)) && all(ismember(expectedDirs, {pastTreeDir.name}))
            %copy only *.gff and roary_output folder
            SourceDir=[pipevar.ref_db_path{1} dl 'Cases' dl 'tree' SS.Case{i} dl];
            CurrDir=[Path.Alignment dl 'tree' SS.Case{i} dl];
            system(['cp ' SourceDir '*.gff ' CurrDir]);
            system(['cp -r ' SourceDir 'roary_output ' CurrDir ]);

        end
    end
end