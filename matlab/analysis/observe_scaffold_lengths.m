if ~exist('Case', 'var')
    load('alignmentReports/all_alignments20_report.mat', 'Case')
end
longestScafs=cellfun(@max, {Case.AssemblyLength});
figure(1);hist(longest)
%there are 13 with counts of 24k or less. what do I do about them?
%maybe remove all scafs of cases that seem contaminated in some way

isPure=cellfun(@length,{Case.GeneName})<5500;
figure(2);hist(longestScafs(isPure))
%no difference
mean(longestScafs(~isPure))
mean(longestScafs(isPure))
%pretty much the same 
%put threshold in 25K to include the majority of scaffolds

%check distribution of ALL scafs
allScafs=cellfun(@(cl) maxk(unique(cl),2), {Case.AssemblyLength}, 'UniformOutput', false)
allScafs=vertcat(allScafs{:});
figure(3);hist(allScafs);


