#!/bin/bash
############################################################################
# ASSEMBLE ALL GENOMES USING SPADES
# Made by: Olga Snitser, 19 May 2020
# Modified by: Noga Aharony, 7 June, 2020
############################################################################
# Paths to input/output files and folders


############################################################################
# Run spades

READSPATH=$PWD/Filtered_data/

my_spades () {

	
	SPADES=/media/kishonylab/KishonyStorage/Apps/SPAdes-3.14.1-Linux/bin/spades.py
	
	if [ ! -f $1/SPAdes_v314_assembly/scaffolds.fasta ]

	then

		rm -rf $1/SPAdes_v314_assembly/
		echo $1

		$SPADES \
		--isolate \
		--threads 1 \
		--pe1-1 $1/R1_trimmed.combined.fastq.gz \
		--pe1-2 $1/R2_trimmed.combined.fastq.gz \
	    -o $1/SPAdes_v314_assembly
	fi
}
export -f my_spades

 # for foldername in $READSPATH/*
 # do
 # 	my_spades "$foldername"
 # done

 parallel --jobs 2 my_spades ::: $READSPATH/*