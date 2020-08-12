#!/bin/bash
############################################################################
# ANNOTATE ALL ASSEMBLED E coli GENOMES USING PROKKA
# Olga Snitser
# 25 May 2020
# Modified by Noga Aharony, June 8, 2020
############################################################################
# Paths to input/output files and folders

#READSPATH=$PWD/Mathews_processing/Filtered_data/


############################################################################
# Run prokka

abricate(){

    BASENAME=$(basename $1)
    echo $1
    echo $2
    if [ ! -f "$2/${BASENAME}/${BASENAME}.gbk" ]
    then

	      echo "abricate ./$1/contigs.fasta --db plasmidfinder > ./$2/$BASENAME.tab" >>abricate_plasmid_batch.sh
	    


	fi



}

export PATH="/media/kishonylab/KishonyStorage/Apps/ncbi-blast-2.10.0+/bin:$PATH"
export -f abricate

read -r READSPATH OUTDIR THREADS <<<$(echo "$1 $2 $3")


parallel --jobs 2 abricate ::: $READSPATH/* ::: $OUTDIR