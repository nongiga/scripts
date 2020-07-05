%% Arrange .gff files for prokka
% Create a roary folder for each 'tree' based on 1000SNP cutoff and copy 
% prokka .gff files from the correct folders into roary 'tree' folders
dl=filesep;
pathT = '/home/kishonystud/kishonyserver/noga/MaccabiUTI/All_same_strain_pairs.xlsx' ;
T = readtable(pathT,'PreserveVariableNames',true) ;
dest_path = '/home/kishonystud/kishonyserver/noga/MaccabiUTI/roary' ;
source_path = '/home/kishonystud/kishonyserver/mathew/' ;

for i = 1:height(T)
   if ~isempty(T.RandomID(i))
       ftree = ['tree' num2str(T.RandomID(i))] ;
      if ~exist([dest_path filesep ftree],'dir')
          mkdir([dest_path filesep ftree])
      end
      pos1=erase(T.Seq_plate_position_1{i},'.');
      pos1=insertBefore(pos1, find(isletter(pos1), 1), '_');
      
      pos2=erase(T.Seq_plate_position_2{i},'.');
      pos2=insertBefore(pos2, find(isletter(pos2), 1), '_');
      
      if exist([source_path dl 'prokka' dl 'Sample_Maccabi_Ecoli_SeqPlate' pos1 '_prokka'], 'dir')
          locname1='Sample_Maccabi_Ecoli_SeqPlate';
      else
           locname1='Sample_Sample_Maccabi_Ecoli_SeqPlate';
      end
      
      
      copyfile([source_path dl 'prokka' dl locname1 pos1 '_prokka' dl 'SeqPlate' pos1 '.gff'],...
                    [dest_path filesep ftree filesep 'Sample_Maccabi_Ecoli_SeqPlate' pos1 '.gff'])
               
      if exist([source_path dl 'prokka' dl 'Sample_Maccabi_Ecoli_SeqPlate' pos2 '_prokka'], 'dir')
          locname2='Sample_Maccabi_Ecoli_SeqPlate';
      else
           locname2='Sample_Sample_Maccabi_Ecoli_SeqPlate';
      end
                
    copyfile([source_path dl 'prokka' dl locname2 pos2 '_prokka' dl 'SeqPlate' pos2 '.gff'],...
            [dest_path filesep ftree filesep 'Sample_Maccabi_Ecoli_SeqPlate' pos2 '.gff'])
   end
end
