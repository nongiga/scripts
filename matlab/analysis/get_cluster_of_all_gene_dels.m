function cllist=get_cluster_of_all_gene_dels(genelist, sg,varargin)
    cllist=[];
    if nargin>2
        curatedlist=varargin{1};
    else
        curatedlist=genelist;
    end
    loc=startsWith(curatedlist(:,1),sg);
    canum=[curatedlist{loc, 11}];
    clnum=[curatedlist{loc, 13}];
    coord=unique([canum' clnum'], 'rows');
    if ~isempty(coord)
    cllist=genelist(ismember(cell2mat(genelist(:,[11 13])),coord, 'rows'),:);
    end
end