function plot_clusters(Path, moptions, pipevar, m, ClusterCase, Case)
    dl=filesep;
    %are there more insertions or deletions?

    %interestingly more deletions

    fig=figure(1);
    for iCase=1:numel(ClusterCase)
        iCase

        Cl=ClusterCase(iCase);
        myCase=Case(Cl.CaseNum);
        del=arrayfun(@(i) Cl.Loc(i,1):Cl.Loc(i,2), 1:size(Cl.Loc,1), 'UniformOutput',false);
        del=[del{:}];
% 
%         genes=vertcat(Cl.Genes{:});
%         ass_num=cellstr(num2str(myCase.AssemblyNum(del)));
%         genes=join([genes, ass_num]);


        [mn, mnidx]=min(myCase.NCov, [],2);
        [mx, mxidx]=max(myCase.NCov, [], 2);    

        subplot(2,1,1)
        hold on
        x=mn';y=mx';
        x(2,:)=nan;
        y(2,:)=nan;

        H = plot(x, y, 'g.');
        arrayfun(@(i) set(H(i), 'Color', [0 0 1]), del);
        %set(H,'ButtonDownFcn',@upon_click)

        %H2 = plot(mn(del), mx(del), 'b.');
        title(sprintf('Num of del = %d/ %d', length(del), numel(mn)))

        subplot(2,1,2)
        hold on
        

        pos=1:size(x,2);
        pos(2,:)=nan;
        H2= plot(pos, x./y, 'g.');

        set(H2,'ButtonDownFcn',@upon_click)
        arrayfun(@(i) set(H2(i), 'Color', [0 0 1]), del);


        set(gca,'Xtick',pos(1,del),'Xticklabel',myCase.GeneName(del),'XtickLabelRotation',90)
        set(gca,'TickLength',[0.001, 0.01])
       pause;
       clf
    end


    %     
    function upon_click(s,~)
        %disp(s.UserData) %
        i=s.SeriesIndex;
        fprintf("%f %f %f %f %f\n",i, mn(i), mx(i), myCase.GeneLength(i), myCase.AssemblyCov(i) )
        %plot_alignment_report(allPhages(m), hMin, hMax, t,ax2);

    end


end