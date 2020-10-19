%% Arrange .gff files for prokka
% Create a roary folder for each 'tree' based on 1000SNP cutoff and copy 
% prokka .gff files from the correct folders into roary 'tree' folders

SS=SameStrains(SameStrains.InstructionsRef==ins,:);
fprintf("copying over prokka files\n");
for i = 1:height(SS)
    %make dir
    sp=SS.seqPlates{i};

    TreeName=['tree' SS.Case{i}];
    TreeFolder=dc_dir(Path.Alignment, TreeName);
    
    
    %if this is a redo and intersection between lists is perfect copy over
    if pipevar.redo(ins)
        
        % set up to see if roary folder was already made (and this is a redo)    
        SourceDir=[pipevar.ref_db_path{ins} dl 'Cases' dl TreeName dl];
        pastTreeDir = dir(SourceDir);
        pastTreeDir(ismember( {pastTreeDir.name}, {'.', '..'}) | [pastTreeDir.isdir]==0) = []; 
        expectedDirs=[strcat(GlobalName, sp), {'roary_output'}];
        
        if isempty(setdiff(expectedDirs, {pastTreeDir.name}))
            %copy only *.gff and roary_output folder

            system(['cp -u ' SourceDir '*.gff ' TreeFolder]);
            system(['cp -r -u ' SourceDir 'roary_output ' TreeFolder ]);
        end
        
    else
        cellfun(@(s) system(['cp -u ' Path.Prokka dl GlobalName s  dl ...
                 GlobalName s  '.gff ' TreeFolder dl ...
                 GlobalName s '.gff']), sp, 'UniformOutput', 0);
    end

end
