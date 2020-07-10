
dl=filesep;
load('gap_data');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';

ins=1;
m=moptions{pipevar.report_multi(ins)+1};
id=[ num2str(pipevar.bp(ins)) m  ];  


Case=struct('Num',{},'GeneName',{},'NCov',{}, 'AssemblyCov', {},'GeneSource', {},'GeneLength',{});

for ins=1:1
    
    id=[ num2str(pipevar.bp(ins)) moptions{pipevar.report_multi(ins)+1} ];

     SS=SameStrains(SameStrains.InstructionsRef==ins,:);
    for i=1:height(SS)
        TreeName=['tree' SS.Case{i}]
        load([Path.Reports dl TreeName '_' id '.mat'], 'myCase')
        for fn = fieldnames(Case)'
           Case(i).(fn{1}) = myCase.(fn{1});
        end
    end
    save([Path.Reports  dl 'all_alignments' id 'trun_report.mat'], 'Case')
end