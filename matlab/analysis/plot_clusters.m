function plot_clusters(Path,pipevar, m, ClusterCase, Case, ins)
    dl=filesep;
    id=[ num2str(pipevar.bp(ins)) m  ];  

    

    fig2=figure(2);
    ax3 = axes(fig2);
    ax3.YLabel.String = 'Normalized Cummulative Sum of Reads';
    ax3.XLabel.String='Location on Gene';
    tl2=sprintf('Gene Name:');ax3.Title.String=tl2;
    axis(ax3,[0 60000 0.01 2])
    grid(ax3)
    hold(ax3)
    
    bc_i=0;
    fig=figure(1);
    
    for iCase=1:numel(ClusterCase)
        
        Cl=ClusterCase(iCase);
        myCase=Case(Cl.CaseNum);
        del=arrayfun(@(i) Cl.Loc(i,1):Cl.Loc(i,2), 1:size(Cl.Loc,1), 'UniformOutput',false);
        del=[del{:}];

        [mn, mnidx]=min(myCase.NCov, [],2);     [mx, mxidx]=max(myCase.NCov, [], 2);    
        x=mn';              y=mx';
        x(2,:)=nan;         y(2,:)=nan;
        pos=1:size(x,2);    pos(2,:)=nan;        

        ax1=subplot(2,1,1);
        H = plot(ax1, x, y, 'g.');
        arrayfun(@(i) set(H(i), 'Color', [0 0 1]), del);
        title(sprintf('Coverage of Each Gene identified in Case # %d', Cl.CaseNum))
        ylabel(ax1, 'Maximum Normalized Coverage')
        xlabel(ax1,'Minimum Normalized Coverage')
        lg=legend(H(del(1)), 'Deleted gene');

        ax2=subplot(2,1,2)
        H2= plot(ax2, pos, x./y, 'g.');
        set(H2,'ButtonDownFcn',@upon_click)
        arrayfun(@(i) set(H2(i), 'Color', [0 0 1]), del);
%         set(ax2,'Xtick',pos(1,:),'Xticklabel',myCase.GeneName,...
%             'XtickLabelRotation',90, 'TickLength',[0.001, 0.01])
        title(sprintf('Ratio of Gene Coverage Along Genome for Case # %d', Cl.CaseNum))
        xlabel(ax2, 'Linearized Gene Number')
        ylabel(ax2, 'Ratio of Maximum and Minimum Coverage')
        pause;
        clf
    end
  
    function upon_click(s,~)
        %disp(s.UserData) %
        i=s.XData(1);
        if ~isfield(myCase, 'GeneLength')
            load(['alignmentReports/tree' num2str(Cl.Num) '_20.mat'], 'myCase')
        end
        fprintf("%f %f %f %f %f\n",i, mn(i), mx(i), myCase.GeneLength(i), myCase.AssemblyCov(i) )

        disp(myCase.GeneName{i});
        
        cNum=Cl.CaseNum;
        gl=myCase.GeneLength(i);
        
        
        cla(ax3)
        for j=1:size(myCase.BaseCov,2)
            plot(ax3,1:gl,cumsum(myCase.BaseCov{i,j})/sum(myCase.BaseCov{i,j}),'-','linewidth',1);
            hold on;            
        end

        
        axis(ax3,[0 myCase.GeneLength(i) 0.01 1])
        legend(ax3, {'Minimum Coverage', 'Maximum Coverage'})
        tl=sprintf("Reads Alignment to Gene %s in Case %6.0f", myCase.GeneName{i}, cNum);
        title(ax3, tl, 'Interpreter', 'none')
        %ax3.Title.String=tl;

    end


end