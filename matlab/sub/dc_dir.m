function mydir=dc_dir(mdir, subdir)
%dc_dir DEFINES and CREATES dir, returns the dirname
    mydir=[mdir filesep subdir];
    if ~exist(mydir, 'dir')
        mkdir(mydir);
    end
end

