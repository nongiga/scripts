dl=filesep;

AD='alignmentReports';
CD='deletion_reports';

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

iCase=2;
    
C=Case(iCase);


[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);

del=mx>0.5 & mn<0.05 ;
sum(del);
sig=C.pval==0;

figure(iCase);

subplot(3,1,1)
hold on
plot(mn, mx, 'g.');
plot(mn(sig), mx(sig), 'r.')
plot(mn(del), mx(del), 'b.')

subplot(3,1,2)
hold on

load([AD dl 'tree' num2str(C.Num) '_20bp'] )

genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')

plot(pos(sig), mn(sig)./mx(sig),'r.')
plot(pos(del), mn(del)./mx(del),'b.')
%set(gca,'Xtick',1:numel(del),'Xticklabel',myCase.GeneName,'XtickLabelRotation',90)
genelength=myCase.GeneLength;
subplot(3,1,3)
hold on
plot(genelength, mn./mx, 'g.')
plot(genelength(del), mn(del)./mx(del),'r.')



%is it gene length that is responsible for all the false positives?
% look at gene length in dels vs no dels
figure;
subplot(2,1,1)
hold on
histogram(genelength, 0:20:1600)
title('gene lengths of all')

subplot(2,1,2)
histogram(genelength(del), 0:20:1600)
title('gene length of deleted genes')

% for C.Num=1001239, there is a very dense cluster of genes w/ mn=0 mx>0.5.
% Did we put the threshold in the right place?

figure;
histogram(mx(del & mn==0), 11)

%plot fold change distribution
figure;
histogram(mn./mx)
hold on
histogram(mn(del)./mx(del), 'FaceColor', 'r')
title('fold change distribution of mn/mn for all and dels')

%the coverage seems to variable in some places
for iCase = 1: numel(Case)

    load([AD dl 'tree' num2str(Case(iCase).Num) '_20bp'] )
    Case(iCase).Reads=myCase.Reads;
    
    
end


x=[];y=[];

%look at the ratio etween minmax over the genome comparison.
if ~exist([AD dl 'mean_vs_pval'])
    figure(1);
    for iCase = 1: 36

        C=Case(iCase);


        [mn, mnidx]=min(C.NCov, [],2);
        [mx, mxidx]=max(C.NCov, [], 2);

        del=mx>0.5 & mn<0.05 ;
        sum(del)


        y=[y nanmean(mn./mx)];
        x=[x sum(C.pval==0)];

        if iCase<37

            load([AD dl 'tree' num2str(C.Num) '_20bp'] )
            sig=myCase.pval==0;
            subplot(6,6,iCase)
            sm=sum(myCase.Reads);
            title([num2str(log10(sm(1))) ' ' num2str(log10(sm(2)))])
            hold on
            genesource=myCase.GeneSource;
            pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
            
            plot(pos, mn./mx, 'g.')
            plot(pos(del), mn(del)./mx(del),'b.')
            plot(pos(sig & del), mn(sig &del )./mx(sig &del),'r.')
            
        end

    end
    save([AD dl 'mean_vs_pval'], 'x', 'y');
    saveas(gcf, ['figures' dl 'minmax_on_genome'])
else
    load([AD dl 'mean_vs_pval'], 'x', 'y')
end


%observe differences in coverage
figure(7);
for iCase = 1: 36

    C=Case(iCase);


    [mn, mnidx]=min(C.NCov, [],2);
    [mx, mxidx]=max(C.NCov, [], 2);

    del=mx>0.5 & mn<0.05 ;
    sum(del)


    y=[y nanmean(mn./mx)];
    x=[x sum(C.pval==0)];

    if iCase<37

        load([AD dl 'tree' num2str(C.Num) '_20bp'] )
        sig=myCase.pval==0;
        subplot(6,6,iCase)
        sm=sum(myCase.Reads);
        title([num2str(log10(sm(1))) ' ' num2str(log10(sm(2)))])
        hold on
        genesource=myCase.GeneSource;
        pos=cellfun(@(x) str2num(x(end-4:end)),genesource);

        plot(mx, mn, 'g.')
        plot(mx(del), mn(del),'b.')
        plot(mx(sig & del), mn(sig &del ),'r.')

    end

end


 %does this issue come from the mn or the mx?
%look at distribution of NCOV values

%it's not clear that one of these is good and one of these is bad.
%let's just take 2 extreme cases and look at their differences
v=arrayfun(@(c) var(var(c.NCov)), Case)
figure (4)

%look at distribution of variance
subplot(2,1,1)
hist(v, 500)
title('distribution of variance')

%is there a connection between variance and the mean of the mn./mx?
subplot(2,1,2)
plot(v,y, '.')
title('variance vs average min/max');

%not necessarily. there are some values with low variance and low
%mean(mn./mx)


%look closer. 
fl=find(y<0.5)
C=Case(fl(1))
figure(6)
subplot(2,1,1)
plot(C.NCov(:,1))
tt1=sprintf("Normalized coverage over each gene in isolate for Case %s", num2str(C.Num));
title(tt1);
hold on
plot(C.NCov(:,2), 'r')
%wow, this is a very extreme example but it really looks like certain large
%parts were deleted

%wait this is one of the 8000-gene ones. time to explore that.
%it looks like the blue has all the genes, so this is not an error in
%classification. but something was added.

%how many genes have an above 0 coverage in each?
sum(C.NCov(:,1)>0.1)
sum(C.NCov(:,2)>0.1)


% in this case it looks like the high variation between min and max is that
% so much is deleted that a lot more reads go to the genome with the
% deleted segments

%what about something more inbetween?
fu=find(y<0.81 & y>0.806)
C=Case(fu(1))


subplot(2,1,2)
plot(C.NCov(:,1))
tt2=sprintf("Normalized coverage over each gene in isolate for Case %s", C.Num);
title(tt2);
hold on
plot(C.NCov(:,2), 'r')


%what is going on with genes in the pangenome that have 0 coverage over all
%the isolates?
zcov=all(myCase.NCov==0,2);
myCase.GeneName(zcov)
myCase.GeneSource(zcov)

%if the reason is due to multialignment each of these genes will get hits
%from the multialigned20.fastq file

 