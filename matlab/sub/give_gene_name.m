function givenNames=give_gene_name(geneName, geneDesc)

% give names based on NAME in descriptions
    isGroup = cellfun(@(s)startsWith(s,'group'),geneName);
    isHypot=contains(geneDesc, 'hypothetical protein');
    nNameyDesc=find(isGroup & ~isHypot);
    lastWord=cellfun(@extractLast, geneDesc(nNameyDesc));
    LWlen=cellfun(@length, lastWord);
    isGeneName=LWlen > 1 & ...
        ((cellfun(@(s) any(regexp(s ,'[0-9]')), lastWord)) | ...
        (LWlen >= 3 & LWlen <=4)) & ...
        ~cellfun(@(s) contains(s, 'ase'), lastWord);
    geneName(nNameyDesc(isGeneName))=lastWord(isGeneName);
    
    
    % clean name underscores and whatnot
    isGroup = cellfun(@(s)startsWith(s,'group'),geneName);
    underscoreLoc = cellfun(@(s) min([find(s=='_',1),inf]), geneName);
    isRemove = underscoreLoc>=4 & underscoreLoc<inf & ~isGroup;
    geneName(isRemove) = cellfun(@(s) extractBefore(s,'_'),geneName(isRemove),'UniformOutput',false);

    %% if same description as a different, named protein, apply it to the unnamed proteins
    [u_desc, ~, ud_idx]=unique(geneDesc);
    for i=1:max(ud_idx)
        sameDesc = ud_idx==i;
        if nnz(sameDesc) > nnz(sameDesc & isGroup)
            newGeneName=unique(lower(geneName(sameDesc & ~isGroup)));
            %issue: the names cannot be applied because there is a number of
            %proteins with same description, different names
            if length(newGeneName)>1
                %if all these names are the same if you remove anything
                %following an underscore...
                isRem=contains(newGeneName,'_');
                clearNewGeneName=newGeneName;
                clearNewGeneName(isRem)=extractBefore(clearNewGeneName(isRem), '_');
                clearNewGeneName=erase(clearNewGeneName, ["orf", "gp"]);
                clearNewGeneName(ismember(clearNewGeneName, 'np'))={'n'};

                sameOneLetter=length(unique(clearNewGeneName))==1;
                %or start in the same 3 letters
                if all(cellfun(@length, newGeneName)>=3)
                    firstThree=cellfun(@(g) g(1:3), newGeneName, 'UniformOutput', 0);
                    startSame=length(unique(firstThree))==1;
                else
                    startSame=0;
                end

                if sameOneLetter || startSame || length(unique(clearNewGeneName))==1
                    geneName(sameDesc & isGroup)= clearNewGeneName(1);
                end

            else
                geneName(sameDesc & isGroup)= newGeneName;
            end

        end
    end

    givenNames=geneName;
    
    function s=extractLast(string)
    sp=split(string, [" ", "-"]);
    s=sp(end);
end
end