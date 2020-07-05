#!/bin/bash
############################################################################
# splits sequnece to 20bp and 
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
			
			if [ ! -f R1_spl20.combined.trimmed.fastq.gz ] && [ -f $SEQPATH/R1_combined.trimmed.fastq.gz ]
			then
				echo "$SEQDIR"

				~/bin/seqkit sliding \
				-s 20 -W 20 \
				$SEQPATH/R1_combined.trimmed.fastq.gz | \
				gzip -c >R1_spl20.combined.trimmed.fastq.gz
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

		echo "$SEQDIR"
		cd $1/$SEQDIR
		# if [ ! -f aligned20.sam ]
		# then 
		# 	bowtie \
		# 	-S -p 1 -v 1 -m 1 \
		# 	--max multialigned20.fastq \
		# 	--un unaligned20.fastq \
		# 	$1/pangenome \
		# 	R1_spl20.combined.trimmed.fastq.gz \
		# 	aligned20.sam #outputq
		# fi

		if [ ! -f aligned_sorted20.bam ]
		then 
			samtools view -bS -o aligned20.bam aligned20.sam
			samtools sort  aligned20.bam -o aligned_sorted20.bam
			samtools index aligned_sorted20.bam
			rm aligned20.bam aligned20.sam
		fi
		
		if [ ! -f aligned20.vcf ]
		then
			samtools mpileup -t AD -ugf $1/roary_output/pan_genome_reference.fa aligned_sorted20.bam > pileup20.vcf.temp
			/media/kishonylab/KishonyStorage/Apps/bcftools/bcftools/bcftools view -o aligned20.vcf pileup20.vcf.temp
			# another line here
			rm pileup20.vcf.temp
		fi

		if [ ! -f depth_clean20.csv ]
		then
			awk -F'[=\t;]' '{print $1 , $2, $9}' aligned20.vcf | sed '/^#/d' > depth_clean20.csv

		fi
	done
}
export -f bowtie_align
echo $READSPATH


echo "split sequence"
parallel --jobs 10 split_sequence ::: $READSPATH/*
echo "bowtie_build"
parallel --jobs 10 bowtie_build ::: $READSPATH/*
echo "bowtie_align"
parallel --jobs 10 bowtie_align ::: $READSPATH/*
