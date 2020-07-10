function cl=get_cluster_location(mn,mx,delet, an, gs)
    cl=[(findstr([0 delet'], [0 1]))' (findstr([delet' 0], [1 0]))'];
    
    %expand clusters if genes around still lowly expressed
    %for each cluster
    len=numel(delet);
    scl=size(cl,1);
    for c=scl:-1:1
        coords=cl(c,:);

        %define range within the genes can be extended:
        %both constrained by beginning of next cluster and range of genes
        if c<scl
            rrange=cl(min(scl,c+1),1)-1;
        else
            rrange=min(len, cl(end)+2);
        end
        
        if c>1
            lrange=cl(max(1,c-1),2)+1;
        else
            lrange=max(1, cl(1)-2);
        end
        generange=lrange:rrange;
        
        % what if there is a 0 coverage seperating 2 clusters?

        rend=intersect([coords(2):coords(2)+4], generange);
        lend=intersect([coords(1)-4:coords(1)], generange);
        
        %further constrain rend and lend to only include genes in the same
        %assembly
        rend=rend(an(rend)==an(coords(2)));
        lend=lend(an(lend)==an(coords(1)));
        
        %diff in mnmx is significant 
        rmnmx=mn(rend)./mx(rend) < 0.5 & mx(rend)>0;
        lmnmx=mn(lend)./mx(lend) < 0.5 & mx(lend)>0;
        
        
        %look for last coordinate in this range that adheres to the current
        %standard
        
        rnewcord=min([find(diff(rmnmx)==-1,1,'first') numel(rmnmx)]);
        
        %from the left: the last increase from 0 to 1 that is followed by 1
        %or if all 1's set as 1 (there will always be at least 1 true)
        
        lnewcord=[find(diff(lmnmx)==1,1,'last')+1 numel(lmnmx)];
        lnewcord=lnewcord(1);
        
        cl(c,2)=rend(rnewcord);
        cl(c,1)=lend(lnewcord);

    end

    
    %loop again to join clusters
    %if clusters are seperated by a single gene (end of 1 +2>start of 2)
    %and are in the same assembly (an(end of 1)==an(start of 2)
    %then join them into one
    temp_cl=[];
     c=1;
     add_last=true;
    while c<scl
        tc=c;
        new_cl=cl(c,:);
        while cl(tc+1,1)<=cl(tc,2)+3 && an(cl(tc+1,1))==an(cl(tc,2))
            tc=tc+1;
            new_cl(2)=cl(tc,2);
            if tc>=scl
                add_last=false;
                break
            end
            
        end
        temp_cl=[temp_cl; new_cl];
        c=tc+1;
    end
    %merge clusters if overlap/1 apart and next one is not in different an
    if add_last
        temp_cl=[temp_cl; cl(c,:)];
    end
    
    cl=temp_cl;
    
    
    %split clusters if placed across different assemblies
    %are the beginning and end in the same assembly?
    dff_ass=find(diff(an(cl)'));
    for da=flip(dff_ass)
        diff_num=diff(an(cl(da,1):cl(da,2)));
        gs_num=unique(gs(cl(da,1):cl(da,2)));
        assert(numel(gs_num)==1, 'error: cluster spread across multiple assemblies')
        %diff_src=diff(gs(cl(da,1):cl(da,2)))
        sploc=find(diff_num);
        spnum=cl(da,1)+sploc-1;
        cl(da+numel(sploc):end+numel(sploc),:)=cl(da:end,:);
        cl(da:da+numel(sploc)-1,2)=spnum;
        cl(da+1:da+numel(sploc),1)=spnum+1;
        %also if assembly is form different source
    end
    c=1;
    while c<= size(cl,1)
        coords=cl(c,:);
        %if cluster has less than 2 genes that have a coverage of zero, remove   
        mnmx2=sort(mn(coords(1):coords(2))./mx(coords(1):coords(2)));
        is_mnmx2=any(mnmx2(1:min(2, length(mnmx2)))>0.0005);
        
        mx2=sort(mx(coords(1):coords(2)));
        ismx2=all(mx2(max(1, length(mx2)-1):end)<0.5);
        
        if  is_mnmx2 || ismx2 || numel(mnmx2)==1
            cl(c,:)=[];
        else
            c=c+1;
        end
    end
        
        
    
    
    
end

