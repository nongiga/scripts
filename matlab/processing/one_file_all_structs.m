load('gap_data.mat');
dl=filesep;
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
        load(FileName, 'myCase')
        if ~exist('Case', 'var')
            Case=myCase;
        end
        Case(i)=myCase;
    end
    
    f=fieldnames(Case)';
    c=cell(length(f),1)
    sCase=cell2struct(c,f);
    for fn = f
       sCase.(fn{1}) = {Case.(fn{1})};
    end
    
    save([Path.Reports  dl 'all_alignment_struct' id '_report.mat'], '-struct','sCase','-v7.3')
end