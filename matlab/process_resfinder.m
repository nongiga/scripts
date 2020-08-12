
ResFinder=readtable('resfinder.tsv', 'FileType', 'text');
load('gap_data', 'SameStrains')
%%


isoName=erase(ResFinder.IsolateID, ["Sample_Maccabi_Ecoli_SeqPlate", "_contigs"]);
plateloc=cellfun(@(in) cellfun(@(sp) ismember(sp, in), ...
     SameStrains.seqPlates, 'Uniformoutput', 0), isoName, 'Uniformoutput', 0);
plateloc=cellfun(@(pl) find(cellfun(@any, pl)), plateloc, 'uniformoutput', 0);

isEmpt=cellfun(@isempty,plateloc);
CaseNum=cell(size(isEmpt));
CaseNum(~isEmpt)=SameStrains.Case(cat(1,plateloc{~isEmpt}));

[toLoad, ia, ic]=unique(CaseNum(~isEmpt));
for i=1:numel(toLoad)

load(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');

fullLoc=find(~isEmpt);
ResDet=ResFinder((fullLoc(ic==i)), [1 2 7 8 9]);
ResDet.Isolate=isoName(fullLoc(ic==i));
ResStruct=table2struct(ResDet(:,[2 4 5 6]), 'toscalar',1);

myCase.Res=ResStruct;

AssemblyData=split(ResDet.Contig, '_',2);
AssemblyData=num2cell(str2double(AssemblyData(:,[2 4 6])),1);
[myCase.Res.Num myCase.Res.Length myCase.Res.Cov]=AssemblyData{:};

save(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');

%for each unique gene in case only one isolate will fit

end