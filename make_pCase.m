function pCase=make_pCase(ClusterCase)
    if ~exist('Case','var')
        load('alignmentReports/all_alignments20_report.mat', 'Case')
    end    
    
    %[Case.pangneome]=NaN; [Case.cGeneIDs]=NaN;
    pCase=struct('pangenome', {}, 'cGeneID',{}, 'Num', {});
    for iCase=1:numel(CLusterCase)
        % pangenome = {'gene', 'description', isPlasmid, isTranspsone, isGroup, hasName} % table of ~4000 genes
        mC=Case(ClusterNum(iCase));

        isPlasmid= num2cell(any(contains(lower(mC.GeneDescription), 'plasmid') | mC.NCov>4, 2));
        isTransposon=num2cell(any(contains(lower(mC.GeneDescription), 'transpos') | startsWith(mC.GeneName, 'tra'), 2));
        isGroup = num2cell(arrayfun(@(s)startsWith(s,'group'),mC.GeneName));
        givenName=give_gene_name(mC.GeneName, mC.GeneDescription);
        hasName=num2cell(~arrayfun(@(s)startsWith(s,'group'),givenName));
        pCase(iCase).Num=iCase;
        pCase(iCase).pangenome=[ mC.GeneName(:), mC.GeneDescription, isPlasmid, isTransposon, isGroup,hasName, givenName, isInsert];

        %if this case contains deletions/insertions
        genelocs=ClusterCase(iCase).Loc;
        pCase(iCase).cGeneID=arrayfun(@(i) [genelocs(i,1):genelocs(i,2)],1:size(genelocs,1), 'UniformOutput',0);

    end
    disp('saving')
    save('all_alignments20_pangenome.mat', 'pCase','-v7.3')
    clear Case
    disp('finished saving')
end