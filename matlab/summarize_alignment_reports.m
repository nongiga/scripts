Case=struct('Num',{},'SeqPlate',{},'GeneName',{},'Date',{},'NCov',{}, 'pval', {}, 'Reads',{});
Path.Reports=dc_dir(Path.Main, 'alignmentReports');

for ins=1:h
    id=[ num2str(pipevar.bp(ins)) moptions{pipevar.report_multi(ins)+1} ];
     SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i=1:height(SS)
        disp(i)
        TreeName=['tree' SS.Case{i}];
        FileName=[Path.Alignment dl TreeName dl 'alignments' id '.mat'];
        %if ~exist([Path.Main dl 'alignmentReports' dl Sub{i} '_20bp.mat'], 'file')
        copyfile(FileName, [Path.Reports dl TreeName '_' id '.mat'])
        load(FileName, 'myCase')
        for fn = fieldnames(Case)'
           Case(i).(fn{1}) = myCase.(fn{1});
        end
    end
    save([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')
end

