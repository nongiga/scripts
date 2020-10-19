
ResFinder=readtable('resfinder.tsv', 'FileType', 'text');
load('gap_data', 'SameStrains')
%%


isoName=erase(ResFinder.IsolateID, ["Sample_Maccabi_Ecoli_SeqPlate", "_contigs", "Sample_"]);
plateloc=cellfun(@(in) cellfun(@(sp) ismember(sp, in), ...
     SameStrains.seqPlates, 'Uniformoutput', 0), isoName, 'Uniformoutput', 0);
plateloc=cellfun(@(pl) find(cellfun(@any, pl)), plateloc, 'uniformoutput', 0);

isEmpt=cellfun(@isempty,plateloc);
CaseNum=cell(size(isEmpt));
CaseNum(~isEmpt)=SameStrains.Case(cat(1,plateloc{~isEmpt}));

[toLoad, ia, ic]=unique(CaseNum(~isEmpt));
fullLoc=find(~isEmpt);
for i=1:numel(toLoad)
    
load(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');

if isfield(myCase, 'Res'), continue, end

disp(toLoad{i})


ResDet=ResFinder((fullLoc(ic==i)), [1 2 3 7 8 9]);
ResDet.Isolate=isoName(fullLoc(ic==i));

AssemblyData=split(ResDet.Contig, '_',2);
AssemblyData=num2cell(str2double(AssemblyData(:,[2 4 6])),1);
[ResDet.Num ResDet.Length ResDet.Cov]=AssemblyData{:};

if iscellstr(ResDet.Start)
    ResDet.Start=str2double(ResDet.Start);
end

if iscellstr(ResDet.End)
    ResDet.End=str2double(ResDet.End);
end  

locs=cellfun(@(sp) ismember(isoName(fullLoc(ic==i)),sp), myCase.SeqPlate, 'UniformOutput', false);

for j=1:numel(myCase.SeqPlate)
    myCase.Res(j)=table2struct(ResDet(locs{j},2:end), 'toscalar',1);
end





save(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');

%for each unique gene in case only one isolate will fit

end