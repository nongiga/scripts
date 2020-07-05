function [baseFilenames, fullFilenames, theFiles] = get_files(myFolder,myPattern)
%GET_FILES Summary of this function goes here
%   Detailed explanation goes here
    filePattern = fullfile(myFolder, myPattern);
    theFiles = dir(filePattern);
    baseFilenames={theFiles.name};
    fullFilenames=arrayfun(@(x) fullfile(myFolder, x), baseFilenames);
    
end

