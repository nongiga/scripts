clear all; close all;
load('gap_data');
Path.Reports='alignmentReports';
Path.Clusters='clusterReports';
dl=filesep;
load([Path.Reports dl 'all_alignments20_report.mat'])
load([Path.Clusters dl 'all_clusters20_report.mat'])
load('cluster_table.mat');
ct=cluster_table;
%%

%MGE length
figure(1);
subplot(2,1,1)
hist(ct.Length, [0:1000:80000])
xlabel('Length of Cluster (bp)')
ylabel('count')


%# genes in MGE
subplot(2,1,2);
hist(ct.NumGenes, 0:1:80);
xlabel('Number of Genes in Cluster')
ylabel('count')

% embedded scaffolds
figure(2);



caseInLength=arrayfun(@(c) sum(c.Length(c.Insert==1)), ClusterCase)
caseDelLength=arrayfun(@(c) sum(c.Length(c.Insert==0)), ClusterCase)
scatter(caseInLength, caseDelLength)
nnz(logical(caseInLength) & ~logical(caseDelLength))
nnz(~logical(caseInLength) & logical(caseDelLength))
nnz(logical(caseInLength) & logical(caseDelLength))
numel(caseDelLength)

meanMinPh=cellfun(@mean, ct.MinPh)
maxMinPh=cellfun(@mean, ct.MaxPh)

%the max loss is of 0.2 in our way to measure average resistance, but gain
%is up to 1.
hist(meanMaxPh(isMeasured)-meanMinPh(isMeasured))

%there are 22 clusters that come with an added resistance above 0.25
nnz(meanMaxPh(isMeasured)-meanMinPh(isMeasured)>0.25)

% think about this measurement of resistance

%only one of the clusters included a gene identified as conferring
%resistance
ct(meanMaxPh(isMeasured)-meanMinPh(isMeasured)>0.25,:).IsRes

%are there difference sin length between genes containing and not
%containing resistance


[h,p]=fishertest([nnz(ct.IsRes & ~ct.Insert) nnz(ct.IsRes & ct.Insert); ...
nnz(~ct.IsRes & ~ct.Insert) nnz(~ct.IsRes & ct.Insert)])

a=histc(ct.Length(ct.IsRes==1), 0:500:15000)
b=histc(ct.Length(ct.IsRes==0), 0:500:15000)
figure
bar(b, 'r')
hold on; bar(a, 'b')
mean(ct.Length(ct.IsRes==1))
mean(ct.Length(ct.IsRes==0))
ttest(ct.Length(ct.IsRes==0), ct.Length(ct.IsRes==1))
ttest(a,b)
ttest2(ct.Length(ct.IsRes==0), ct.Length(ct.IsRes==1))
[h,p]=ttest2(ct.Length(ct.IsRes==0), ct.Length(ct.IsRes==1))