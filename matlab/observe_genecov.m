close all
clear all

dl=filesep;
load('path')
AD=[Path.Main dl 'alignmentReports'];

load([AD  dl 'all_alignments20bp_report.mat'], 'Case')

CD=[Path.Main dl 'deletion_reports'];
figure(1);


x=[];y=[];

%look at the ratio etween minmax over the genome comparison.
for iCase = 1: numel(Case)
    
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
        title([ strjoin(C.SeqPlate) nanmean(mn./mx)])
        hold on
        genesource=myCase.GeneSource;
        pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
        plot(pos, mn./mx, 'g.')
        plot(pos(sig), mn(sig)./mx(sig),'r.')
        plot(pos(del), mn(del)./mx(del),'b.')
    end
    
end
save([AD dl 'mean_vs_pval'], 'x', 'y');
%saveas(gcf, [Path.Main dl 'figures' dl 'minmax_on_genome'])


%in some figures there is a huge variation between minmax and in some there
%isn't. can we tie it back to the overall coverage or the variance of
%coverage over a genome?

figure; 
%looks like when plotting # of sigs vs mean of mn./mx, lower mean gives more sigs
subplot(2,1,1)
plot(x,y, '.')

%when plotting y there are 2 distributions
subplot(2,1,2)
hist(y,500)
 
%at the hist there is a dip at around ~0.81. Is that a clear border between
% the 'good' and the 'bad'-quality  mn./mx graphs?

%this is empty
sum(y>0.81 & y<0.811)

%there are a few locs here
sum(y>0.811 & y<0.8114)
fa=find(y>0.811 & y<0.8114)
C=Case(fa(1));

fu=find(y<0.81 & y>0.806)

figure(3)
%here, at an average mn./mx of ~0.811 it seems quite clear

load(['tree' num2str(C.Num) '_20bp'] )
[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);
del=mx>0.5 & mn<0.05 ;
sig=myCase.pval==0;
subplot(2,1,1)
title([ strjoin(C.SeqPlate) nanmean(mn./mx)])
hold on
genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')
plot(pos(sig), mn(sig)./mx(sig),'r.')
plot(pos(del), mn(del)./mx(del),'b.')


%what about a little under?
C=Case(fu(1));
subplot(2,1,2)
load(['tree' num2str(C.Num) '_20bp'] )
[mn, mnidx]=min(C.NCov, [],2);
[mx, mxidx]=max(C.NCov, [], 2);
del=mx>0.5 & mn<0.05 ;
sig=myCase.pval==0;
figure(3)
title([ strjoin(C.SeqPlate) nanmean(mn./mx)])
hold on
genesource=myCase.GeneSource;
pos=cellfun(@(x) str2num(x(end-4:end)),genesource);
plot(pos, mn./mx, 'g.')
plot(pos(sig), mn(sig)./mx(sig),'r.')
plot(pos(del), mn(del)./mx(del),'b.')
%does this issue come from the mn or the mx?
%look at distribution of NCOV values

%it's not clear that one of these is good and one of these is bad.
%let's just take 2 extreme cases and look at their differences
v=arrayfun(@(c) var(var(c.NCov)), Case)
figure (4)

%look at distribution of variance
subplot(2,1,1)
hist(v, 500)

%is there a connection between variance and the mean of the mn./mx?
subplot(2,1,2)
plot(v,y, '.')

%not necessarily. there are some values with low variance and low
%mean(mn./mx)


%look closer. 
fl=find(y<0.5)
C=Case(fl(1))
figure(5)
subplot(2,1,1)
plot(C.NCov(:,1))
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
C=Case(fu(1))
subplot(2,1,2)
plot(C.NCov(:,1))
hold on
plot(C.NCov(:,2), 'r')


%what is going on with genes in the pangenome that have 0 coverage over all
%the isolates?
zcov=all(myCase.NCov==0,2);
myCase.GeneName(zcov)
myCase.GeneSource(zcov)

%if the reason is due to multialignment each of these genes will get hits
%from the multialigned20.fastq file

