#!/bin/bash
############################################################################
# CREATE PANGENOMES FOR ALL 18 'TREES' (1000SNP CUTOFF) USING ROARY
# Olga Snitser
# 04 June 2020
# Modified by Noga Aharony, June 8, 2020
############################################################################
# Paths to input/output files and folders
READSPATH=$PWD/roary

############################################################################
# Run roary

split_sequence(){
		for SP in $1/*.gff ; do
			SEQDIR=$(basename -- $SP)
			SEQDIR="${SEQDIR%.*}"

			SEQPATH=~/MaccabiUTI/Mathews_processing/Filtered_data/$SEQDIR/
			mkdir -p $1/$SEQDIR; cd $1/$SEQDIR
			
			if [ ! -f R1_spl40.combined.trimmed.fastq.gz ] && [ -f $SEQPATH/R1_combined.trimmed.fastq.gz ]
			then
				echo "$SEQDIR"

				~/bin/seqkit sliding \
				-s 60 -W 60 -g \
				$SEQPATH/R1_combined.trimmed.fastq.gz | \
				~/bin/seqkit seq -m 30 | \
				gzip -c >R1_spl40.combined.trimmed.fastq.gz
			fi
		done

}

export -f split_sequence

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




bowtie_align() {
	for SP in $1/*.gff ; do
		SEQDIR=$(basename -- $SP)
		SEQDIR="${SEQDIR%.*}"


		cd $1/$SEQDIR
		if [ ! -s aligned.sam ]
		then 
			bowtie \
			-S -p 1 -v 3 -m 1 \
			--max multialigned.fastq \
			--un unaligned.fastq \
			$1/pangenome \
			R1_spl40.combined.trimmed.fastq.gz \
			aligned.sam #outputq
		fi

		if [ ! -s hits.ls ]
		then
			awk '{print $3}' aligned.sam | \
			 sort | \
			 uniq -c | \
			 sort -nr >hits.ls 
		fi
	done
}

export -f bowtie_align
echo $READSPATH



parallel --jobs 6 split_sequence ::: $READSPATH/*
parallel --jobs 6 bowtie_build ::: $READSPATH/*
parallel --jobs 15 bowtie_align ::: $READSPATH/*
