function create_case_alignment(IsoNum,seqPlates,RawSubDirName,CaseNum,SampleDate, Path, GlobalName, bp) 

dl=filesep;

SPLoc=cellfun(@(t) find(endsWith(RawSubDirName, sprintf('SeqPlate%s', t))),seqPlates, 'UniformOutput', false);
SPLoc=[SPLoc{:}];
nSPs=length(seqPlates);


%get pangenome info into myCase.Gene
Fasta_Info=fastaread([Path.Alignment dl 'tree' num2str(CaseNum) dl 'roary_output' dl 'pan_genome_reference.fa']);
sp=regexp({Fasta_Info.Header}','\s','split','once');
nGenes=length(sp);

%get gene origin location 


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
    'Assembly', {cell(nGenes, 1)},...
    'Strand', {repmat('.',nGenes, 1)},...
    'Reads', {zeros(nGenes, nSPs)});

%initialize base coverage matrix
myCase.BaseCov=repmat(arrayfun(@(gl) zeros(gl,1), myCase.GeneLength, 'UniformOutput', false),1,nSPs);
myCase.NormalizedCov=myCase.BaseCov;
clear Fasta_Info sp

assert(issorted(myCase.GeneSource),'error: Source names not sorted. Can cause issues later.');

for s=1:nSPs

    IsolateDir=[Path.Alignment dl 'tree' num2str(myCase.Num) dl GlobalName myCase.SeqPlate{s}]
    
    %get prokka detail
    %get info from .gff file to get location of strand on assembly
    AtGFF=GFFAnnotation([IsolateDir '.gff']);
    At=getData(AtGFF, [1:find(cellfun(@isempty, AtGFF.Attributes), 1, 'first')-1]);
    
    as=cellfun(@(at) extractBetween(at(1:25), '=',';'), {At.Attributes}, 'UniformOutput', false);
    as=[as{:}];   
    %assign data from .gff struct to the Case struct via merging tables
    
    %I can only do this idx assignment since both lists are already sorted
    [~, ia,ib]=intersect(myCase.GeneSource, as);
    myCase.GeneStart(ia)=[At(ib).Start];
    myCase.GeneEnd(ia)=[At(ib).Stop];
    Assembly=split({At(ib).Reference}, '_');
    myCase.AssemblyNum(ia)=str2double(Assembly(1, :, 2));
    myCase.AssemblyLength=str2double(Assembly(1,:, 4));
    myCase.AssemblyCov=str2double(Assembly(1,:, 6));

    myCase.Strand(ia)=At(ib).Strand;
    
    fid=fopen([IsolateDir dl 'depth_clean'  num2str(bp) '.csv']);
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
%normalized for all genes in an isolate (find median basecov on a gene and
%then median these)

%reads per gene
myCase.Reads=cellfun(@(bc) sum(bc), myCase.BaseCov)/100;
%sum of reads
SR=sum(myCase.Reads);

%calculate normalized cov expected
k=(myCase.GeneLength-bp)*SR./(sum(myCase.GeneLength)-bp*numel(myCase.GeneLength));
myCase.NCov=myCase.Reads./k;

%min and max of rate of reads gene takes
[mn, mnidx]=min(myCase.NCov, [],2);
[mx, mxidx]=max(myCase.NCov, [], 2);

%the number of reads in min and max
mnR=arrayfun(@(i) myCase.Reads(i, mnidx(i)), 1:length(mnidx))';
mxR=arrayfun(@(i) myCase.Reads(i, mxidx(i)), 1:length(mxidx))';

%calculate z score, pval, stderr
stderr=sqrt((mnR./(SR(mnidx)'.^2))+(mxR./SR(mxidx)'.^2));
err=(mx-mn);
z=err./stderr;
myCase.pval=normcdf(z, 'upper');

m='';
if pipevar.multia

save([Path.Alignment dl 'tree' num2str(myCase.Num) dl 'alignments'  bp 'bp'], 'myCase');

