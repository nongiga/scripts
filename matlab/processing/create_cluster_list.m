Path.Clusters='clusterReports';
dl=filesep;
load([Path.Clusters  dl 'all_clusters20_report.mat'], 'ClusterCase')
cluster_table=[];

for iCase = 1: numel(ClusterCase)
    C=ClusterCase(iCase);
    C.CaseNum=repmat(C.CaseNum, size(C.Genes));
    C.Num=repmat(str2double(C.Num), size(C.Genes));
    C.NumGenes=C.NumGenes';
    C.Description=C.Description';
    C.Start=C.Loc(:,1)';
    C.End=C.Loc(:,2)';

    
    C.MinPh=num2cell(zeros(7,numel(C.MinDate)),1);
    isResMin=any(datenum(C.PhDate)==datenum(C.MinDate),1);
    
    C.MinPh(isResMin)=arrayfun(@(i) cell2mat(struct2cell(C.PhMes(C.PhDate==i'))), C.MinDate(isResMin), 'UniformOutput',0);
    
    
    isResMax=any(datenum(C.PhDate)==datenum(C.MaxDate),1);
    C.MaxPh=num2cell(zeros(7,numel(C.MaxDate)),1);
    
    %Trimethoprim/Sulfa, Ciprofloxacin, Amoxicillin/CA, Cefuroxime - Axetil, Cephalexin, Nitrofurantoin, Fosfomycin
    C.MaxPh(isResMax)=arrayfun(@(i) cell2mat(struct2cell(C.PhMes(C.PhDate==i'))), C.MaxDate(isResMax), 'UniformOutput',0);
    
    
    C=structfun(@transpose, C, 'UniformOutput',0);

    
    C=rmfield(C,{'N', 'Loc', 'PhMes','PhPlate','PhDate'});
    cluster_table=[cluster_table; struct2table(C)];
    

end

cluster_table

cluster_table.Genes=cellfun(@transpose, cluster_table.Genes, 'UniformOutput', false);
cluster_table.Description=cellfun(@transpose, cluster_table.Description, 'UniformOutput', false);
cluster_table.ResGenes=cellfun(@transpose, cluster_table.ResGenes, 'UniformOutput', false);
cluster_table.ResPh=cellfun(@transpose, cluster_table.ResPh, 'UniformOutput', false);

isMeasured=~any([cluster_table.MinPh{:}; cluster_table.MaxPh{:}]==0)';
resDiff=cellfun(@minus, cluster_table.MaxPh(isMeasured), cluster_table.MinPh(isMeasured), 'UniformOutput', false);
%isMeasured=find(isMeasured)
%does resistance increase upon insertion?
sum([resDiff{:}],2)


save('cluster_table','cluster_table')
ct=cluster_table;
isType= ct.IsPhage | ct.IsPlasmid | ct.IsTra;
ct(~isType& ct.IsRes,:).Genes{:};

%are insertions enriched for resistance genes?
[h,p]=fishertest([nnz(ct.IsRes & ~ct.Insert) nnz(ct.IsRes & ct.Insert); ...
            nnz(~ct.IsRes & ~ct.Insert) nnz(~ct.IsRes & ct.Insert)])


% what to do nanopore on

%cases where phenotypic resistance increased upon insertion
resDiff=cellfun(@minus, cluster_table.MaxPh, cluster_table.MinPh, 'UniformOutput', false);
resIncrease=any([resDiff{:}]>1,1)';


%cases that had identified resistance genes
[ct(resIncrease & isMeasured & ct.Insert & ct.IsRes & ~ct.IsPlasmid,:).MaxPh{:}]

%does an insertion of a resistance gene correlate with presence of
%resistance?
a=nnz(ct.IsRes & isMeasured & ~resIncrease)
b=nnz(ct.IsRes & isMeasured & resIncrease)
c=nnz(~ct.IsRes & isMeasured & ~resIncrease)
d=nnz(~ct.IsRes & isMeasured & resIncrease)

% yes
[h, p]=fishertest([a b ; c d])

%does it get down to the resultion of fitting antibiotic to its measured
%resisatance?