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

READSPATH=$PWD/Filtered_data/

my_atropos() {
	#trim adapters
	
	FWDCONT=CTGTCTCTTATA
	REVCONT=$FWDCONT
	EXT=combined.fastq
	ATROPOS=/media/kishonylab/KishonyStorage/Apps/anaconda3/bin/atropos

	
	if [ ! -f $1/R1_trimmed.$EXT ] && [ ! -f $1/R2_trimmed.$EXT ] && [ -f $1/R1_$EXT ] && [ -f $1/R2_$EXT ]
	
	then
		echo $1/R1_$EXT

		$ATROPOS trim \
		--threads 4 \
		-pe1 $1/R1_$EXT \
		-pe2 $1/R2_$EXT \
		-a $FWDCONT \
		-A $REVCONT \
		--error-rate 0.1 \
		--insert-match-error-rate 0.2 \
		--quality-cutoff 15,15 \
		--minimum-length 50 \
		--max-n 1 \
		-o $1/R1_trimmed.$EXT \
		-p $1/R2_trimmed.$EXT
	fi

}

 export -f my_atropos
 

 # for foldername in $READSPATH/*
 # do
 # 	my_atropos "$foldername" 
 # done


 parallel --jobs 20 my_atropos ::: $READSPATH/*