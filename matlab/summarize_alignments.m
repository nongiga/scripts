function summarize_alignments
dl=filesep;
load('gap_data', 'Path')

Case=struct('Num',{},'SeqPlate',{},'GeneName',{},'Date',{},'NCov',{}, 'pval', {});
AD=[Path.Main dl 'alignmentReports'];

%if ~exist([AD  dl 'all_alignments20bp_report.mat'], 'file')
i=1;
while 1<=length(Sub)
    disp(i)
    FileName=[Path.Alignment dl Sub{i} dl 'alignments20bp.mat'];
    %if ~exist([Path.Main dl 'alignmentReports' dl Sub{i} '_20bp.mat'], 'file')
    if exist(FileName, 'file')
        copyfile(FileName, [Path.Main dl 'alignmentReports' dl Sub{i} '_20bp.mat'])
        load(FileName, 'myCase')
        for fn = fieldnames(Case)'
           Case(i).(fn{1}) = myCase.(fn{1});
        end
        i=i+1;
    else
        Sub(i)=[];
    end
    %end


end

    save([AD  dl 'all_alignments20bp_report.mat'], 'Case')
    save(['matfiles' dl 'path'], 'Path', 'Sub')
%end








