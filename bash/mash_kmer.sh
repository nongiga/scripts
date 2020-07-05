#!/bin/bash
############################################################################
# MASH_KMER ANALYSIS
# Noga Aharony, June 8, 2020
# Based on Mathew Stracy's scripts
############################################################################
# trims all files to appropriate length for Kmer analysis
############################################################################


# count number of reads in file to determine lengh to truncate
DATA_PATH=$PWD/Filtered_data/
FOLDER_KEY="Sample_Maccabi_Ecoli_SeqPlate"


for FOLDERNAME in $DATA_PATH/*
do
echo $(cat $FOLDERNAME/*trimmed* |wc -l)/4|bc
done


mash_kmer(){
	FOLDERNAME=$1
	#sample a subset of the reads. Mathew used 500000 so I guess it is sufficient
	seqtk sample -s11  $FOLDERNAME/R1_trimmed.combined.fastq 500000 > $FOLDERNAME/R1_trimmed_500K.fastq
	#create kmer map
    mash sketch -r -m 2 -s 500000 $FOLDERNAME/R1_trimmed_500K.fastq
}


export -f mash_kmer

parallel --jobs 10 mash_kmer ::: $DATA_PATH/*

#paste all sequences together
KMER_DIR=$PWD/mash_kmers_analysis
mkdir KMER_DIR
mash paste $KMER_DIR/all_isolates Filtered_data/$FOLDER_KEY*/*.fastq.msh
mash dist all_isolates.msh all_isolates.msh > mash_dist_all_isolates_500K