load('seq_res')
load('clusterReports/all_clusters20_report.mat')
Path.Clusters='clusterReports';
DrugNames={'Trimethoprim', 'Ciprofloxacin', 'Amoxicillin', 'Cefuroxime', 'Cephalexin', 'Nitrofurantoin', 'Fosfomycin'};

%Trimethoprim/Sulfa, Ciprofloxacin, Amoxicillin/CA, Cefuroxime - Axetil, Cephalexin, Nitrofurantoin, Fosfomycin
seq_res.Properties.VariableNames{'DrishaNum'} = 'IsoNum';
 [u, uidx]=unique(IsolatesNames.IsoNum);
CleanIsolatesNames=(IsolatesNames(ismember(1:height(IsolatesNames), uidx)' & ~isnan(IsolatesNames.IsoNum),:));
IsolatesProfiles=innerjoin(seq_res,CleanIsolatesNames, 'Keys',{'IsoNum'});

%add numerical resistance to cluster case
locs=str2double({ClusterCase.Num})==IsolatesProfiles.NewRandomID;

IsolatesProfiles=[IsolatesProfiles array2table(IsolatesProfiles.numerical_res, ...
    'VariableNames',DrugNames)];

clusterRes=arrayfun(@(i) {table2struct(IsolatesProfiles(locs(:,i),DrugNames))}, 1:size(locs,2));
clusterPlate=arrayfun(@(i) {erase(IsolatesProfiles.SeqSubDirName(locs(:,i),:),'Sample_Maccabi_Ecoli_SeqPlate')}, 1:size(locs,2));
clusterDate=arrayfun(@(i) {IsolatesProfiles.SampleDate(locs(:,i),:)}, 1:size(locs,2));


[ClusterCase.PhMes]=clusterRes{:};
[ClusterCase.PhPlate]=clusterPlate{:};
[ClusterCase.PhDate]=clusterDate{:};
save([Path.Clusters  dl 'all_clusters20_report.mat'], 'ClusterCase')


% clusterDates