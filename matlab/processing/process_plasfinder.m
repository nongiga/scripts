
PlasFinder=readtable('plasfinder.tsv', 'FileType', 'text');
load('gap_data', 'SameStrains')

%%


pisoName=erase(PlasFinder.IsolateID, ["Sample_Maccabi_Ecoli_SeqPlate", "_contigs", "Sample_"]);
pplateloc=cellfun(@(in) cellfun(@(sp) ismember(sp, in), ...
     SameStrains.seqPlates, 'Uniformoutput', 0), pisoName, 'Uniformoutput', 0);
pplateloc=cellfun(@(pl) find(cellfun(@any, pl)), pplateloc, 'uniformoutput', 0);

pisEmpt=cellfun(@isempty,pplateloc);
pCaseNum=cell(size(pisEmpt));
pCaseNum(~pisEmpt)=SameStrains.Case(cat(1,pplateloc{~pisEmpt}));

[toLoad, ia, ic]=unique(pCaseNum(~pisEmpt));
for i=1:numel(toLoad)
    
    load(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');
    disp(toLoad{i})


    fullLoc=find(~pisEmpt);
    PlasDet=PlasFinder((fullLoc(ic==i)), [1 2 6 7 8 9]);
    PlasDet.Isolate=pisoName(fullLoc(ic==i));
    PlasStruct=table2struct(PlasDet(:,[2:6]), 'toscalar',1);
    myCase.Plas=PlasStruct;
    if iscellstr(myCase.Plas.Start)
        myCase.Plas.Start=str2double(myCase.Plas.Start);
    end
    if iscellstr(myCase.Plas.End)
        myCase.Plas.End=str2double(myCase.Plas.End);
    end  
    
    AssemblyData=split(PlasDet.Contig, '_',2);
    AssemblyData=num2cell(str2double(AssemblyData(:,[2 4 6])),1);
    [myCase.Plas.Num myCase.Plas.Length myCase.Plas.Cov]=AssemblyData{:};

    save(['alignmentReports/tree' toLoad{i} '_20.mat'],'myCase');

    %for each unique gene in case only one isolate will fit

end