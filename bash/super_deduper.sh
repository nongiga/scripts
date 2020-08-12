#!/bin/bash
############################################################################
# ATROPOS
# Made by Olga Snitser, 19 May 2020
# Modified by Noga Aharony, June 7, 2020
############################################################################
# Paths to input/output files and folders


############################################################################
# adapter contamination to remove

############################################################################
# Run Atropos

READSPATH=/media/kishonylab/KishonyStorage/noga/MaccabiUTI/Mathews_processing/raw_seqs/

my_superdeduper() {
	#trim adapters


	NODUPPATH=/media/kishonylab/KishonyStorage/noga/MaccabiUTI/super_deduper/
	BASENAME=$(basename -- $1)
	BASENAME=$BASENAME.nodup

	if [ ! -f $NODUPPATH/$BASENAME.nodup_R1.fastq.gz ]; then
		if [ -f $1/*_L001_R1_001.fastq.gz ]; then
			hts_SuperDeduper \
			-1 $1/*_L001_R1_001.fastq.gz \
			-2 $1/*_L001_R2_001.fastq.gz \
			-f $NODUPPATH/$BASENAME
		elif [ -f $1/*_L006_R1_001.fastq.gz ]; then
			hts_SuperDeduper \
			-1 $1/*_L006_R1_001.fastq.gz \
			-2 $1/*_L006_R2_001.fastq.gz \
			-f $NODUPPATH/$BASENAME
		elif [ -f $1/*_L008_R1_001.fastq.gz ]; then
			hts_SuperDeduper \
			-1 $1/*_L008_R1_001.fastq.gz \
			-2 $1/*_L008_R2_001.fastq.gz \
			-f $NODUPPATH/$BASENAME
		fi
	fi

}

 export -f my_superdeduper
 

 # for foldername in $READSPATH/*
 # do
 # 	my_atropos "$foldername" 
 # done


 parallel --jobs 3 my_superdeduper ::: $READSPATH/*
