function nl = count_lines(fn)

nls = evalc (['!wc -l ' fn]) ;
ft = strfind(nls, ' ') ;
nl = str2double(nls(1:ft-1)) ;