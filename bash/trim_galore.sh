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

NODUPPATH=/media/kishonylab/KishonyStorage/noga/MaccabiUTI/super_deduper/

my_trimgalore() {
	#trim adapters

	OUTDIR=/media/kishonylab/KishonyStorage/noga/MaccabiUTI/trim_galore/
	NODUPPATH=/media/kishonylab/KishonyStorage/noga/MaccabiUTI/super_deduper/
	BASENAME=$(basename -- $1)
	BASENAME="${BASENAME%.*.*.*}"
	echo $OUTDIR/$BASENAME
	if [ ! -d $OUTDIR/$BASENAME ]; then
			trim_galore --fastqc --paired  \
			$NODUPPATH/$BASENAME.nodup_R1.fastq.gz \
			$NODUPPATH/$BASENAME.nodup_R2.fastq.gz \
			-o $OUTDIR/$BASENAME/
	fi

}

 export -f my_trimgalore
 

 # for foldername in $READSPATH/*
 # do
 # 	my_atropos "$foldername" 
 # done


parallel --jobs 10 my_trimgalore ::: $NODUPPATH/*
