function pCase=make_pCase(ClusterCase, Case)
   
    phage_genes={'bet','Alpha','S_','N_','C_','4_','3_','int','cox',...
        'B_','Q_','O_','L_','X_','Y_','R_','V_','ant','Alpha','ninB',...
        'P_','cII','gam','bet','exo', 'yokD'};

    ClusterNum=[ClusterCase.CaseNum];
    %[Case.pangneome]=NaN; [Case.cGeneIDs]=NaN;
    pCase=struct('pangenome', {}, 'cGeneID',{}, 'Num', {});
    for iCase=1:numel(ClusterCase)
        % pangenome = {'gene', 'description', isPlasmid, isTranspsone, isGroup, hasName} % table of ~4000 genes
        mC=Case(ClusterNum(iCase));
        
        isPlasmid= num2cell(any(contains(lower(mC.GeneDescription), 'plasmid') | mC.NCov>4, 2));
        plasmidScaffold=unique(mC.AssemblyNum([isPlasmid{:}]));
        isPlasmidCorrect=num2cell(ismember(mC.AssemblyNum,plasmidScaffold));

        isTransposon=num2cell(any(contains(lower(mC.GeneDescription), 'transpos') | startsWith(mC.GeneName, 'tra'), 2));
        
        isGenomic=num2cell(mC.AssemblyLength>120000);
        
        isPhage=num2cell(startsWith(mC.GeneName, phage_genes));
        
        isGroup = num2cell(arrayfun(@(s)startsWith(s,'group'),mC.GeneName));
        givenName=give_gene_name(mC.GeneName, mC.GeneDescription);
        hasName=num2cell(~arrayfun(@(s)startsWith(s,'group'),givenName));
        
        pCase(iCase).Num=iCase;
        pCase(iCase).pangenome=[ mC.GeneName(:), mC.GeneDescription, ...
            isPlasmid, isTransposon, isGroup,hasName, givenName, ...
            isPlasmidCorrect, isGenomic,isPhage];

        %if this case contains deletions/insertions
        genelocs=ClusterCase(iCase).Loc;
        pCase(iCase).cGeneID=arrayfun(@(i) [genelocs(i,1):genelocs(i,2)],1:size(genelocs,1), 'UniformOutput',0);
        
        % is it an insert
        cInsert=arrayfun(@(i) repmat(ClusterCase(iCase).Insert(i),1,(genelocs(i,2)-genelocs(i,1)+1)),1:size(genelocs,1), 'UniformOutput',0);
        cInsert=[cInsert{:}];
        pCase(iCase).cInsert=cInsert;
        
        % cluster number
        cNum=arrayfun(@(i) repmat(i,1,(genelocs(i,2)-genelocs(i,1)+1)),1:size(genelocs,1), 'UniformOutput',0);
        cNum=[cNum{:}];
        pCase(iCase).cNum=cNum;
        
        %is it a phage        
        cPhage=arrayfun(@(i) repmat(ClusterCase(iCase).IsPhage(i),1,(genelocs(i,2)-genelocs(i,1)+1)),1:size(genelocs,1), 'UniformOutput',0);
        cPhage=[cPhage{:}];
        pCase(iCase).cPhage=cPhage;
        
        
        
    end
    disp('saving')
    save('all_alignments20_pangenome.mat', 'pCase','-v7.3')
    disp('finished saving')
end