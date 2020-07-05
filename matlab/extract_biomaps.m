function extract_biomaps(Path)

dl=filesep;
if ~exist('all_alignments.mat', 'file')
    allAlign=[];
    for c=1:length(Path.Sub.Case)
        if exist([Path.Alignment dl Path.Sub.Case{c} dl 'alignments.mat'],'file')
            load([Path.Alignment dl Path.Sub.Case{c} dl 'alignments.mat']);
        else
            continue
        end

        myPhages= [x.Template.Phage];
        myGenBanks=[myPhages.GenBank];
        allAlign=[allAlign myPhages(([myGenBanks.Category]< 3) & ([myPhages.Length] > 10000))]; %get all phages with category 1/2
    end

    myAlign=[(arrayfun(@(mp) mp.Align(mp.Min),allAlign))' (arrayfun(@(mp) mp.Align(mp.Max),allAlign))'];

    %sort so that earlier date come sfirst
%         
    Dates=reshape([myAlign.Date],[],2);
    toSwitch=Dates(1,:)<Dates(2,:);
    myAlign(toSwitch,:)=[myAlign(toSwitch,2) myAlign(toSwitch,1)];

    %myRate=reshape([myAlign.Count]./[myAlign.NSeqs],[],2);
    myRate=reshape([myAlign.NormalizedCount],[],2);

    save([Path.Main dl 'all_phages'], 'allPhages','myAlign','myRate','-v7.3')
end
    
end