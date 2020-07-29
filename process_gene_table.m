%% process geneTable to give more group names and descriptions

% give names based on descriptions
isHypot=contains(geneTable(:,2), 'hypothetical protein');
nNameyDesc=find(isGroup & ~isHypot)
lastWord=cellfun(@extractLast, geneTable(nNameyDesc,2));
LWlen=cellfun(@length, lastWord);
isGeneName=LWlen > 1 & ...
    ((cellfun(@(s) any(regexp(s ,'[0-9]')), lastWord)) | ...
    (LWlen >= 3 & LWlen <=4)) & ...
    ~cellfun(@(s) contains(s, 'ase'), lastWord);
geneTable(nNyDesc(isGeneName), 1)=lastWord(isGeneName);

% give description based on name
nameSearch=geneTable(~isGroup & isHypot,1);
isRemove=contains(nameSearch,'_');
nameSearch(isRemove)=extractBefore(nameSearch(isRemove), '_');
descSearch=cellfun(@(ns) find(startsWith(geneTable(:,1),ns) & ~isHypot) , nameSearch,'UniformOutput', false)
descLocs=cellfun(@(ds) ds(1), descSearch);
geneTable(~isGroup & isHypot,2)=geneTable(descLocs, 2)

%% group all proteins based on cluster

ClusterT=readtable('blast/undefined_p.fa.groups', 'FileType', 'text');

for i=1:height(ClusterT)
    gene_loc=ismember(entry_name, ClusterT{i,:});
    
    c= ismember(geneTable(gene_loc,2), 'hypothetical protein');
    % if all the proteins are hypothetical change desc
    if all(c)
        clustered_descriptions(gene_loc)={['cluster ' num2str(i)]};
    %if one or more of the genes is defined but the rest are not
    elseif sum(c)
        % apply the description to all hypothetical proteins
        %divide descriptions to 
        descriptions=geneTable(gene_loc & ~isHypot,2);
        disp(descriptions)
        clustered_descriptions(gene_loc & isHypot)=descriptions(1);
        
    end
    
end


function s=extractLast(string)
    sp=split(string, [" ", "-"]);
    s=sp(end);
end