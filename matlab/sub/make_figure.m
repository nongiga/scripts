function make_figure(okgenes, delgenes, sigGenes, iClusterCase)
    cloc=arrayfun(@(sg) horzcat(okgenes{startsWith(okgenes(:,1),sg),11}), sigGenes, 'uniformoutput', 0);
    sigMat=zeros(numel(sigGenes), numel(iClusterCase));
    for i=1:numel(cloc)
        sigMat(i, cloc{i})=1;
    end
 %   gdelgenes=delgenes;
%    delgenes=vertcat(delgenes{:});
    sigDesc=arrayfun(@(sg) okgenes(find(startsWith(okgenes(:,1), sg), 1, 'first'),2), sigGenes);
    % look at the 

    IndelSC=[sum(sigMat,2)];
    figure(2);  h=my_imagesc(IndelSC);
    set(gca, 'ytick',[1:numel(sigGenes)],'yticklabel', join([sigGenes sigDesc], ','))
    set(gca, 'xtick',[1 2],'xticklabel', {'Insertion', 'Deletion'})
    colorbar

    for ii = 1:size(IndelSC,1)
        for jj = 1:size(IndelSC,2)
            set(h(ii,jj),'ButtonDownFcn',@buttondown,'UserData',[ii,jj])
        end
    end

    function buttondown(s, ~)
        ij = s.UserData;
        sg=sigGenes(ij(1));
        cllist=get_cluster_of_all_gene_dels(delgenes, sg,okgenes)
        
    end

end
