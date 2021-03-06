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
conda activate virsorter
export PERL5LIB=/media/kishonylab/KishonyStorage/noga/perl5/lib/perl5/
virsorter(){

	SPADESPATH="spades_assembly"
    BASENAME=$(basename $1)
    echo $BASENAME

    if [ ! -f virsorter/${BASENAME}/Predicted_viral_sequences/VIRSorter_prophages_cat-6.fasta ]
    then
    		rm -rf virsorter/${BASENAME}/
		   ~/opt/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl \
		   -f assembly/${BASENAME}/contigs.fasta \
		   --wdir virsorter/${BASENAME} \
		   --db 1 --ncpu 1 \
		   --data-dir ~/opt/virsorter-data/ 
	fi



}
export -f virsorter

#read -r READSPATH OUTDIR THREADS <<<$(echo "$1 $2 ")


parallel --jobs 6 virsorter ::: assembly/* 
