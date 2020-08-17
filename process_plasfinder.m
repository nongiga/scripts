
PlasFinder=readtable('plasfinder.tsv', 'FileType', 'text');
load('gap_data', 'SameStrains')

%%


isoName=erase(PlasFinder.IsolateID, ["Sample_Maccabi_Ecoli_SeqPlate", "_contigs"]);
plateloc=cellfun(@(in) cellfun(@(sp) ismember(sp, in), ...
     SameStrains.seqPlates, 'Uniformoutput', 0), isoName, 'Uniformoutput', 0);
plateloc=cellfun(@(pl) find(cellfun(@any, pl)), plateloc, 'uniformoutput', 0);

isEmpt=cellfun(@isempty,plateloc);
CaseNum=cell(size(isEmpt));
CaseNum(~isEmpt)=SameStrains.Case(cat(1,plateloc{~isEmpt}));

[toLoad, ia, ic]=unique(CaseNum(~isEmpt));
for i=1:numel(toLoad)
    
    if isfield(myCase, 'Plas'), continue, end
    disp(toLoad{i})


    fullLoc=find(~isEmpt);
    PlasDet=PlasFinder((fullLoc(ic==i)), [1 2 6 7 8 9]);
    PlasDet.Isolate=isoName(fullLoc(ic==i));
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