function create_case_alignment(IsoNum,seqPlates,RawSubDirName,CaseNum,SampleDate, Path, GlobalName, bp, m,r) 

dl=filesep;
id=[ num2str(bp) m ];
nSPs=length(seqPlates);
if r
    %get details of seqplates in isolates list
    SPLoc=cellfun(@(t) find(endsWith(RawSubDirName, sprintf('SeqPlate%s', t))),seqPlates, 'UniformOutput', false);
    SPLoc=[SPLoc{:}];
    

    %get pangenome info into myCase.Gene
    Fasta_Info=fastaread([Path.Alignment dl 'tree' CaseNum dl 'roary_output' dl 'pan_genome_reference.fa']);
    sp=regexp({Fasta_Info.Header}','\s','split','once');
    nGenes=length(sp);

    %initialize my case
    myCase=struct('Num', CaseNum, ...
        'SeqPlate', {seqPlates'}, ...
        'IsNum', {IsoNum(SPLoc)'}, ...
        'Sequence', {{Fasta_Info.Sequence}'}, ...
        'GeneLength', {cellfun(@(gs) length(gs), {Fasta_Info.Sequence}')}, ...
        'GeneSource', {cellfun(@(s) s(1),sp)}, ...
        'GeneName', {cellfun(@(s) s(2),sp)},...
        'BaseCov', {cell(nGenes,nSPs)},...
        'Date', {arrayfun(@(i) {SampleDate(i)}, SPLoc)},...
        'NCov', {zeros(nGenes, nSPs)}, ...
        'GeneStart', {zeros(nGenes, 1)}, ...
        'GeneEnd', {zeros(nGenes, 1)}, ...
        'Strand', {repmat('.',nGenes, 1)},...
        'AssemblyNum', {zeros(nGenes, 1)}, ...
        'AssemblyLength', {zeros(nGenes, 1)}, ...
        'AssemblyCov', {zeros(nGenes, 1)}, ...
        'Reads', {zeros(nGenes, nSPs)},...
        'GeneDescription', {cellfun(@(s) s(1),sp)});

    %initialize base coverage matrix
    myCase.BaseCov=repmat(arrayfun(@(gl) zeros(gl,1), myCase.GeneLength, 'UniformOutput', false),1,nSPs);
    myCase.NormalizedCov=myCase.BaseCov;
 
    for s=1:nSPs

        IsolateDir=[Path.Alignment dl 'tree' num2str(myCase.Num) dl GlobalName myCase.SeqPlate{s}];
        disp(IsolateDir);

        %get prokka detail
        %% get info from .gff file to get location of strand on assembly
        AtGFF=GFFAnnotation([IsolateDir '.gff']);
        At=getData(AtGFF, [1:find(cellfun(@isempty, AtGFF.Attributes), 1, 'first')-1]);

        as=cellfun(@(at) extractBetween(at(1:25), '=',';'), {At.Attributes}, 'UniformOutput', false);
        description=cellfun(@(at) extractAfter(at, 'product='), {At.Attributes}, 'UniformOutput', false)';
        as=[as{:}];   
        %assign data from .gff struct to the Case struct via merging tables

        %I can only do this idx assignment since both lists are already sorted
        [~, ia,ib]=intersect(myCase.GeneSource, as);
        
        myCase.GeneDescription(ia)=description(ib);
        myCase.GeneStart(ia)=[At(ib).Start];
        myCase.GeneEnd(ia)=[At(ib).Stop];
        Assembly=split({At(ib).Reference}, '_');
        myCase.AssemblyNum(ia)=str2double(Assembly(1, :, 2));
        myCase.AssemblyLength(ia)=str2double(Assembly(1,:, 4));
        myCase.AssemblyCov(ia)=str2double(Assembly(1,:, 6));
        myCase.Strand(ia)=At(ib).Strand;

        %% get alignment details
        fid=fopen([IsolateDir dl 'depth_clean' id '.csv']);
        HitsT=textscan(fid, "%s %d %d");
        %for each line in table place value in right location
        g=1;
        for i=1:length(HitsT{1})
            while ~isequal(HitsT{1}(i), myCase.GeneSource(g))
                g=g+1;
            end
            myCase.BaseCov{g,s}(HitsT{2}(i))=HitsT{3}(i);
        end
        fclose(fid);
    end

else
    load([Path.Alignment dl 'tree' CaseNum{:} dl 'alignments' id '.mat'], 'myCase')
    for s=1:nSPs
        IsolateDir=[Path.Alignment dl 'tree' num2str(myCase.Num) dl GlobalName myCase.SeqPlate{s}];
        %get prokka detail
        %% get info from .gff file to get location of strand on assembly
        AtGFF=GFFAnnotation([IsolateDir '.gff']);
        At=getData(AtGFF, 1:find(cellfun(@isempty, AtGFF.Attributes), 1, 'first')-1);

        as=cellfun(@(at) extractBetween(at(1:25), '=',';'), {At.Attributes}, 'UniformOutput', false);
        as=[as{:}];   
        [~, ia,ib]=intersect(myCase.GeneSource, as);
        
        description=cellfun(@(at) extractAfter(at, 'product='), {At.Attributes}, 'UniformOutput', false);
        myCase.GeneDescription(ia)=description(ib);
    end    
end

%reads per gene
myCase.Reads=cellfun(@(bc) sum(bc), myCase.BaseCov)./bp;
k=(myCase.GeneLength-bp)*sum(myCase.Reads)./(sum(myCase.GeneLength)-20*numel(myCase.GeneLength));
myCase.NCov=myCase.Reads./k;

Reads=cellfun(@(bc) sum(bc), myCase.BaseCov)./100;
SR=sum(Reads);
Rate=Reads./SR;

[mn, mnidx]=min(Rate, [],2);
[mx, mxidx]=max(Rate, [], 2);

stderr=sqrt((mn./(SR(mnidx)'))+(mx./SR(mxidx)'));
err=(mx-mn);
z=err./stderr;
myCase.pval=normcdf(z, 'upper');

save([Path.Alignment dl 'tree' myCase.Num dl 'alignments'  id  '.mat'], 'myCase');

