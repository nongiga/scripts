load('gap_data.mat');
dl=filesep;
Case=struct('Num',{},'SeqPlate',{},'GeneName',{},'Date',{},'NCov',{}, 'pval', {}, 'Reads',{});

for ins=1:h
    
    id=[ num2str(pipevar.bp(ins)) moptions{pipevar.report_multi(ins)+1} ];
%     if exist([Path.Reports  dl 'all_alignments' id '_report.mat'], 'file')
%         continue
%     end
    
    
     SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i=1:height(SS)
        disp(i)
        TreeName=['tree' SS.Case{i}];
        FileName=[Path.Alignment dl TreeName dl 'alignments' id '.mat'];
        copyfile(FileName, [Path.Reports dl TreeName '_' id '.mat'])
        load(FileName, 'myCase')
        for fn = fieldnames(Case)'
           Case(i).(fn{1}) = myCase.(fn{1});
        end
    end
    save([Path.Reports  dl 'all_alignments' id '_report.mat'], 'Case')
end

