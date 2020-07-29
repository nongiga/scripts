#!/bin/bash
############################################################################
# CREATE PANGENOMES FOR ALL 18 'TREES' (1000SNP CUTOFF) USING ROARY
# Olga Snitser
# 04 June 2020
# Modified by Noga Aharony, June 8, 2020
############################################################################
# Paths to input/output files and folders
############################################################################
# Run roary
echo "running roary in dir $1 in $2 parallel jobs"

roary_func() {
	FN=$1
	if [ ! -f $FN/roary_output/pan_genome_reference.fa ]
	then
		rm -rf $FN/roary_output/
		echo "$FN"
		roary \
		-o $(basename $FN) \
		-f $FN/roary_output \
		-e --mafft \
		$FN/*.gff
	fi
}

export -f roary_func
echo $READSPATH
parallel --jobs 15 roary_func ::: $1/* ::: $2

