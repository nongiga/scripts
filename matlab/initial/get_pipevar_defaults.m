function defaults = get_pipevar_defaults

%**
if ~exist('pipevar_defaults.xlsx', 'file')
    pathToDefaults = '/media/kishonylab/KishonyStorage/Sequencing_analysis' ; 
else
    pathToDefaults='';
end
% if running locally or not is server, manually change line below to
% reflect where you keep 'pipevar_defaults.xlsx'. For example, if the
% server is mounted locally at /home/user/kishonyserver/' then:
% pathToDefaults = '/home/user/kishonyserver/Sequencing_analysis' ; 
d = readtable([pathToDefaults filesep 'pipevar_defaults.xlsx']) ;
defaults=d;

return