
%CREATE_ALIGNMENT_BIOMAPS Summary of this function goes here
%   Detailed explanation goes here

%   move these out when possible

dl=filesep;

disp("create alignment biomaps")
for ins=1:height(pipevar)
    
    
    m=moptions{pipevar.report_multi(ins)+1};
    id=[ num2str(pipevar.bp(ins)) m ];
    
    mm=pipevar.mm(ins);
    %get same strain isolates and casenums
    IS=IsolatesNames(IsolatesNames.InstructionsRef==ins,:);
    SS=SameStrains(SameStrains.InstructionsRef==ins,:);

    r=cellfun(@(c) ~exist([Path.Alignment dl 'tree' c dl 'alignments' id '.mat'],'file'), SS.Case);
    l=SS.Case; s=SS.seqPlates;
    
    %if it is not a redo, only start making reports from scratch to the
    %isolates for which files do not exist (in a redo will also go over
    %reports and 'fix' them)
    
    if ~pipevar.recreate_reports(ins)
        l=l(r); s=s(r); 
    end
        
    %loop through each case
    if pipevar.parallel(ins)
        parfor (f=1:mm)
            for i=f:mm:numel(l)
                disp(i);
                create_case_alignment(IS.IsoNum,[s{i}],IS.RawSubDirName,l(i),IS.SampleDate, Path, GlobalName, pipevar.bp(ins), m, r(i)) 

            end

        end
    else
        for i=1:1:numel(l)
            disp(i);
            create_case_alignment(IS.IsoNum,[s{i}],IS.RawSubDirName,l(i),IS.SampleDate, Path, GlobalName, pipevar.bp(ins), m, r(i)) 

        end
    end
end

