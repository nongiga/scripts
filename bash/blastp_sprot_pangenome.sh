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

#!/bin/bash
############################################################################
# splits sequnece to $2bp. Requires Seqkit.
#called by gap (matlab function)
# used for MaccabiUTI project
# Noga Aharony
# July 10, 2020
############################################################################
# $1: main dir, $2: bp length, $4: parallel, $3=input sequence folder
export BLASTDB=/media/kishonylab/KishonyStorage/DBs/
blastp_pangenome(){

		cd $1
		echo $1
		if [ ! -f pan_genome_reference_blast_sprot.tsv ] 
		then
			blastp -outfmt "7 qseqid sseqid evalue " -query pan_genome_reference.faa -db uniprot_sprot -out pan_genome_reference_blast_sprot.tsv
		fi

}

export -f blastp_pangenome

echo "blastp pangenome"
read -r maindir threads <<<$(echo "$1 $2")


parallel --jobs $threads blastp_pangenome ::: $maindir/* 
