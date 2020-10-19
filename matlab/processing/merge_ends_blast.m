T=readtable('all_ends.tsv', 'Filetype','text', 'delimiter', '\t');


isSame=arrayfun(@(i) isequal(T.Var1(i),T.Var2(i)), 1:height(T));
T=T(~isSame,:);

clusters={};

for i=1:height(T)
    isIn=0;
    l=T{i,[1 2]};
    for j=length(clusters):-1:1
        if any(ismember(clusters{j},l))
            clusters(j)={unique([clusters{j}  l])};
            isIn=1;
            break;
        end
    end
    if ~isIn
        clusters=[clusters; {l}];
    end
end

B=readtable('Blast_all_ends.html', 'Filetype','text', 'delimiter', '\t');
large=find(cellfun(@numel, clusters)>100)

for l=large'
    iLines=B{ismember(B{:,1},clusters{l}),2};

    a = unique(iLines);

    [s,~,j]=unique(iLines);
    s{mode(j)}
end

