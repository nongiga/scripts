function results = my_blast(seq, query, dbase)


fid = fopen(query,'W') ;
fprintf(fid,'>query1\n%s\n',seq ) ;
fclose(fid) ;
try
    %warning('off','bioinfo:blastreadlocal:NoHits');
      results = blastlocal(['-i ' query ' -d "' dbase '" -p blastn -e 1e-7 -F F -W 32']) ;
catch
    results = [] ;
end
return