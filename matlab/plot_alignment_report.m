function g=plot_alignment_report(myCase, i, hMin, hMax, ax)
%DISPLAY ALIGNMENT REPORT Summary of this function goes here
%   Detailed explanation goes here

    %load alignment file
    seqp=split(seqp);
   %TO DO
   %display title at correct location
   %make overlap clearer
   %display prophage #, name
   %buttons for all 77 pairs
    
   
%    
    alignfolder=[mainfolder 'alignment_reports/alignReport_Maccabi_Ecoli_SeqPlate'];
     if exist([alignfolder seqp{1} '.mat'], 'file') && exist([alignfolder seqp{2} '.mat'], 'file')
         Min=load([alignfolder seqp{1} '.mat'], 'Alignments');
         Min=cell2mat(Min.Alignments{j,1});
         Max=load([alignfolder seqp{2} '.mat'], 'Alignments');
         Max=cell2mat(Max.Alignments{j,1});
         disp(proname);
         %load([sourcefolder '/' proname '.mat'], 'g');
         set(hMin,'XData',1:length(Min),'YData',cumsum(Min)/sum(Min));
         set(hMax,'XData',1:length(Max),'YData',cumsum(Max)/sum(Max));
         axis(ax,[0 length(Max) 0.01 1])
         g=get_ref_genome(proname);
%          cellfun(@(x) disp(x(1,:)), {g.CDS.product}')
%          contains({g.CDS.product}','hypothetical')
%          cds_no_hyp=~contains([g.CDS.product], 'hypothetical');
%          g_wo_hyp=g.C
%          disp()
        tl=sprintf("Reads Alignment to Phage %s on isolate %6.0f", proname, cNum);
        ax.Title.String=tl;

     end
end