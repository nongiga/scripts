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

prokka(){


    BASENAME=$(basename $1)

    if [ ! -f "$2/${BASENAME}/${BASENAME}.gbk" ]
    then

	    if [ ! -f "$2/${BASENAME}_Viruses/${BASENAME}_Viruses.gbk" ]
	    then
		    /media/kishonylab/KishonyStorage/Apps/prokka/bin/prokka \
		    --force \
		    --outdir "$2/${BASENAME}_Viruses" \
		    --prefix "${BASENAME}_Viruses" \
		    --locustag $BASENAME \
		    --evalue 0.05 \
		    --coverage 1 \
		    --cpus 1 \
		    --kingdom Viruses \
		    $1/contigs.fasta
		fi


	      /media/kishonylab/KishonyStorage/Apps/prokka/bin/prokka \
	    --force \
	    --outdir $2/$BASENAME \
	    --prefix $BASENAME \
	    --locustag $BASENAME \
	    --evalue 0.05 \
	    --cpus 1 \
	    --coverage 1 \
	    --kingdom Bacteria \
	    --proteins "$2/${BASENAME}_Viruses/${BASENAME}_Viruses.gbk" \
	    $1/contigs.fasta


	fi



}

export PATH="/media/kishonylab/KishonyStorage/Apps/ncbi-blast-2.10.0+/bin:$PATH"
export -f prokka
echo "running prokka"
read -r READSPATH OUTDIR THREADS <<<$(echo "$1 $2 $3")


parallel --jobs $THREADS prokka ::: $READSPATH/* ::: $OUTDIR