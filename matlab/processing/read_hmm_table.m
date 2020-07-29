clear all; close all;

fid=fopen('pan_genome_hmm.tsv');

varnames={ 'target','t_accession','query' ,'q_accession','b_Evalue',...
    'b_score','b_bias','n_Evalue','n_score','n_bias','exp','reg',...
    'clu','ov','env','dom','rep','inc','description'};

C = textscan(fid,'%s', 'delimiter', '\n');
lines=C{:}(4:end-10);
%skip first 3 lines and last 9 lines
sp=cellfun(@(s) split(s, ' '), lines, 'UniformOutput', 0);
sp=cellfun(@(s) s(~cellfun(@isempty, s)), sp, 'UniformOutput', 0);
sp=cellfun(@(s) [s(1:18)' join(s(19:end))], sp, 'UniformOutput', 0);
T=cell2table(vertcat(sp{:}),'VariableNames',varnames);

load('alignments20.mat')

pfam={'','',''};
for i=1:numel(myCase.GeneSource)
    idx=find(ismember(T.query,myCase.GeneSource(i)));
    if isempty(idx) 
        pfam(i,:)={'','',''};
    else
        [mn, midx]=min(str2double(T.b_Evalue(idx)));
        pfam(i,:)=T{idx(midx), {'target', 't_accession','description'}};
    end
    % target, accession, description of target
end



