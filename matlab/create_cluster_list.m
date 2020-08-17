%load([Path.Clusters  dl 'all_clusters20_report.mat'], 'ClusterCase')
cluster_table=[];

for iCase = 1: numel(ClusterCase)
    C=ClusterCase(iCase);
    C.CaseNum=repmat(C.CaseNum, size(C.Genes));
    C.Num=repmat(str2double(C.Num), size(C.Genes));
    C.GeneNum=C.GeneNum';
    C.Description=C.Description';
    C.Start=C.Loc(:,1)';
    C.End=C.Loc(:,2)';

    
    C.MinPh=num2cell(zeros(7,numel(C.MinDate)),1);
    isResMin=any(datenum(C.PhDate)==datenum(C.MinDate),1);
    
    C.MinPh(isResMin)=arrayfun(@(i) cell2mat(struct2cell(C.PhMes(C.PhDate==i'))), C.MinDate(isResMin), 'UniformOutput',0)
    
    
    isResMax=any(datenum(C.PhDate)==datenum(C.MaxDate),1);
    C.MaxPh=num2cell(zeros(7,numel(C.MaxDate)),1);
    
    C.MaxPh(isResMax)=arrayfun(@(i) cell2mat(struct2cell(C.PhMes(C.PhDate==i'))), C.MaxDate(isResMax), 'UniformOutput',0)
    
    
    C=structfun(@transpose, C, 'UniformOutput',0);
    
    C=rmfield(C,{'N', 'Loc', 'PhMes','PhPlate','PhDate'});
    cluster_table=[cluster_table; struct2table(C)];
    

end

cluster_table
save('cluster_table','cluster_table')
ct=cluster_table;
isType= ct.IsPhage | ct.IsPlasmid | ct.IsTra;
ct(~isType& ct.IsRes,:).Genes{:}


