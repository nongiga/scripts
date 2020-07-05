#!/bin/bash
############################################################################
# ANNOTATE ALL ASSEMBLED E coli GENOMES USING PROKKA
# Olga Snitser
# 25 May 2020
# Modified by Noga Aharony, June 8, 2020
############################################################################
# Paths to input/output files and folders

READSPATH=$PWD/Filtered_data/
SPADESPATH="SPAdes_v314_assembly"

############################################################################
# Run prokka

prokka(){
    BASENAME=$(basename $1)
    export PATH="/media/kishonylab/KishonyStorage/Apps/ncbi-blast-2.10.0+/bin:$PATH"
    /media/kishonylab/KishonyStorage/Apps/prokka/bin/prokka \
    --force \
    --outdir $1/prokka \
    --prefix $BASENAME \
     \
    $1/$SPADESPATH/scaffolds.fasta
}


export -f prokka
parallel --jobs 10 prokka ::: $READSPATH/*