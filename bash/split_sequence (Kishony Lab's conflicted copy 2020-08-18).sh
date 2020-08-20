#!/bin/bash
############################################################################
# splits sequnece to $2bp. Requires Seqkit.
#called by gap (matlab function)
# used for MaccabiUTI project
# Noga Aharony
# July 10, 2020
############################################################################
# $1: main dir, $2: bp length, $4: parallel, $3=input sequence folder

split_sequence(){
		for SP in $1/*.gff ; do
			SEQDIR=$(basename -- $SP)
			SEQDIR="${SEQDIR%.*}"
			SEQPATH=$3/$SEQDIR/
			mkdir -p $1/$SEQDIR; cd $1/$SEQDIR
			echo $SEQPATH
			if [[ ! -f R1_spl$2.combined.trimmed.fastq.gz &&  -f $SEQPATH/R1_combined.trimmed.fastq.gz ]]
			then
				seqkit sliding \
				-s $2 -W $2 \
				$SEQPATH/R1_combined.trimmed.fastq.gz | \
				gzip -c >R1_spl$2.combined.trimmed.fastq.gz
			fi
		done

}

export -f split_sequence

echo "split sequence"
read -r maindir inputdir bp threads <<<$(echo "$1 $2 $3 $4")


parallel --jobs $threads split_sequence ::: $maindir/* ::: $bp ::: $inputdir