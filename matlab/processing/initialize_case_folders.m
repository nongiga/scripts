%% Arrange .gff files for prokka
% Create a roary folder for each 'tree' based on 1000SNP cutoff and copy 
% prokka .gff files from the correct folders into roary 'tree' folders

SS=SameStrains(SameStrains.InstructionsRef==ins,:);
for i = 1:height(SS)
    %make dir
    sp=SS.seqPlates{i};

    TreeName=['tree' SS.Case{i}];
    TreeFolder=dc_dir(Path.Alignment, TreeName);
    SourceDir=[pipevar.ref_db_path{ins} dl 'Cases' dl TreeName dl];

    pastTreeDir = dir(SourceDir);
    pastTreeDir(ismember( {pastTreeDir.name}, {'.', '..'}) | [pastTreeDir.isdir]==0) = []; 
    expectedDirs=[strcat(GlobalName, sp), {'roary_output'}];

    %if this is a redo and intersection between lists is perfect copy over
    if pipevar.redo(ins) && isempty(setdiff(expectedDirs, {pastTreeDir.name}))
        %copy only *.gff and roary_output folder

        system(['cp -u ' SourceDir '*.gff ' TreeFolder]);
%             system(['cp -u ' SourceDir '*.faa ' TreeFolder]);
        system(['cp -r -u ' SourceDir 'roary_output ' TreeFolder ]);
    else
        fprintf("%s has changed\n", TreeName);
        for j=1:numel(sp)
            system(['cp -u ' Path.Prokka dl '*' GlobalName sp{j}  '_prokka' dl ...
                'SeqPlate' sp{j} '.gff ' TreeFolder dl GlobalName sp{j} '.gff']);
%                 system(['cp -u ' Path.Prokka dl '*' GlobalName sp{j}  '_prokka' dl ...
%                     'SeqPlate' sp{j} '.faa ' TreeFolder dl GlobalName sp{j} '.faa']);
         end
    end

end
