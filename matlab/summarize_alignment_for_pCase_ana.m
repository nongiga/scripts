load('gap_data.mat');
dl=filesep;
Case=struct('Num',{},'SeqPlate',{},'GeneName',{},'GeneDescription', {},'Date',{},'NCov',{},'Reads',{}, 'AssemblyNum',{}, 'AssemblyLength',{});
Path.Reports='alignmentReports';
Path.Alignment=Path.Reports;

for ins=1:h
    
    id=[ num2str(20) moptions{pipevar.report_multi(ins)+1} ];
     SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i=1:height(SS)
        disp(i)
        TreeName=['tree' SS.Case{i}];
        FileName=[Path.Alignment dl TreeName '_20.mat'];
%         copyfile(FileName, [Path.Reports dl TreeName '_' id '.mat'])
        load(FileName, 'myCase')
        for fn = fieldnames(Case)'
           Case(i).(fn{1}) = myCase.(fn{1});
        end
    end
    save([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')
end