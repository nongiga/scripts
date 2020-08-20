load('cluster_table');

ct=cluster_table
%how much is there of each?
size(ct,1) %=186
nnz(ct.IsPhage) %=73
nnz(ct.IsPlasmid) %=72
nnz(ct.IsTra) %=41

nnz(ct.IsPlasmid & ct.IsPhage) %=27 plasmid & phage?!
nnz(ct.IsPlasmid & ct.IsTra) %=16 plasmid & transposon
nnz(ct.IsPhage & ct.IsTra) %=20 transposable phage

nnz(ct.IsPlasmid & ct.IsPhage & ct.IsTra) %=9

nnz(~ct.IsPlasmid & ~ct.IsPhage & ~ct.IsTra) %=54 none of the above

%within resistance, is there a bias to gain or lose?
nnz(ct.IsRes & ~ct.Insert)
nnz(ct.IsRes & ct.Insert)