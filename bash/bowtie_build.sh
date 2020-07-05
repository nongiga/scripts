#!/bin/bash
############################################################################
# Build database for bowtie
#called by gap (matlab function)
# used for MaccabiUTI project
# Noga Aharony
# July 10, 2020
############################################################################
# $1: maindir $2: # of jobs


bowtie_build(){
	cd $1
	if [ ! -f pangenome.1.ebwt ]
	then
		bowtie-build \
		roary_output/pan_genome_reference.fa \
		pangenome
	fi

}

export -f bowtie_build


echo "bowtie_build"
parallel --jobs $2 bowtie_build ::: $1/*
