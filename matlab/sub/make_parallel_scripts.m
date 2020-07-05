function [subscripts, subscript_names]=make_parallel_scripts(ScriptPath, ScriptName, mm, varargin)

    dl=filesep;
    
   
    %delete to ensure no interference
    delete([ScriptPath dl ScriptName '*.sh']);
    
    %make subscripts & handles
    subscript_names=arrayfun(@(i) [ScriptPath dl ScriptName '_' num2str(i) '.sh'],1:mm,'UniformOutput',false);
    subscripts=arrayfun(@(n) fopen(n{:},'w'), subscript_names);
    arrayfun(@(n) fprintf(n, '#!/bin/bash\n'), subscripts);
    
    %if need be add something to the top of the scripts
    
    if size(varargin,1)==1
        arrayfun(@(n) fprintf(n, [varargin{1} '\n']), subscripts);
    end
    
    %make mainscript and write orders in it
    mainscript_name=[ScriptPath dl ScriptName '.sh'];
    mainscript=fopen(mainscript_name,'w');
    fprintf(mainscript, '#!/bin/bash\n');
    arrayfun(@(s) fprintf(mainscript,'%s & \n',s{:}), subscript_names);
    fclose(mainscript);
    
    %make files executable
    arrayfun(@(sn) fileattrib(sn{:}, '+x'), subscript_names);
    fileattrib(mainscript_name, '+x');
    
    
end