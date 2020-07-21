dl=filesep;
InstructionLines = get_instructions;

h=height(InstructionLines);
%get pipevar defaults
pipevar = get_pipevar_defaults; 
pipevar= repmat(pipevar, h,1);
fn = intersect(InstructionLines.Properties.VariableNames, pipevar.Properties.VariableNames) ;
pipevar(:,fn) = InstructionLines(:,fn) ;

% form path
PathName=pwd;
Path=struct('Main', PathName, ...
            'Alignment',[PathName dl 'Cases'], ...
            'Reports', [PathName dl 'alignmentReports'], ...
            'Clusters', [PathName dl 'clusterReports'],...
             'Prokka', [PathName dl 'prokka'],...
             'Assembly', [PathName dl 'assembly']);
structfun(@mymkdir, Path);

% Getting isolates details and cretaing list of isolates to select from
IsolatesNames = arrayfun(@(i) get_isolates_details(InstructionLines.isolates_file_name{i},i), 1:h,'UniformOutput', false) ;
IsolatesNames=vertcat(IsolatesNames{:});

%unify the sample date and time
SampleDate=arrayfun(@(i) datetime(IsolatesNames.SampleDate(i)), 1:height(IsolatesNames))';
IsolatesNames.SampleDate=SampleDate;

%get global name from files
GlobalName=IsolatesNames.RawSubDirName{1};
GlobalName=GlobalName(1:(regexp(GlobalName, '\d','once')-1));

%set parallel processing to be off if I am on the mounted path from the lab
%computer
pipevar_option=[pipevar.mm(1),0];
pipevar.parallel=repmat(pipevar_option(isfolder('/home/kishonystud/kishonyserver/')+1), h, 1);

% get same strains
%SameStrains=get_same_strains(pipevar.same_strain_isolates, IsolatesNames, GlobalName);
SameStrains=arrayfun(@(i) get_same_strains(InstructionLines.same_strain_isolates{i},IsolatesNames, GlobalName,i), 1:h,'UniformOutput', false) ;
SameStrains=vertcat(SameStrains{:});

moptions={'','.multi'};

% save all collected variables
save([Path.Main dl 'gap_data'], 'Path', 'pipevar', 'IsolatesNames', 'GlobalName', 'SameStrains', 'h', 'moptions');

copy_over_assemblies

function mymkdir(d)
    if ~isfolder(d), mkdir(d); end
end


