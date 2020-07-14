#!/bin/bash
############################################################################
# splits sequnece to $2bp. Requires Seqkit.
#called by gap (matlab function)
# used for MaccabiUTI project
# Noga Aharony
# July 10, 2020
############################################################################
# $1: main dir, $2: bp length, $4: parallel, $3=input sequence folder
hmmsearch_pangneome(){

		cd $1/roary_output
		echo $1
		if [ ! -f pan_gneome_hmm.tsv ] 
		then
			hmmscan --cpu 1 --tblout pan_gneome_hmm.tsv /media/kishonylab/KishonyStorage/DBs/Pfam-A.hmm pan_genome_reference.fa
		fi

}

export -f hmmsearch_pangneome

echo "hmmsearch pangenome"
read -r maindir threads <<<$(echo "$1 $2")


parallel --jobs $threads hmmsearch_pangneome ::: $maindir/* 