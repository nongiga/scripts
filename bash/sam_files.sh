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


bowtie_align() {
	for SP in $1/*.gff ; do
		SEQDIR=$(basename -- $SP)
		SEQDIR="${SEQDIR%.*}"


		cd $1/$SEQDIR
		if [ ! -f aligned.sam ]
		then 
			bowtie \
			-S -p 1 -v 3 -m 1 \
			--max multialigned.fastq \
			--un unaligned.fastq \
			$1/pangenome \
			R1_spl40.combined.trimmed.fastq.gz \
			aligned.sam #outputq
		fi

		if [ ! -f aligned_sorted.bam ]
		then 
			samtools view -bS -o aligned.bam aligned.sam
			samtools sort  aligned.bam -o aligned_sorted.bam
			samtools index aligned_sorted.bam
			rm aligned.bam
		fi
		
		if [ ! -f aligned.vcf ]
		then
			samtools mpileup -t AD -ugf $1/roary_output/pan_genome_reference.fa aligned_sorted.bam > pileup.vcf.temp
			/media/kishonylab/KishonyStorage/Apps/bcftools/bcftools/bcftools view -o aligned.vcf pileup.vcf.temp
			rm pileup.vcf.temp
		fi

		if [ ! -f depth.csv ]
		then
			awk -F'[=\t;]' '{print $1 , $2, $9}' aligned.vcf >depth.csv

		fi

		if [ ! -f depth_clean.csv ]
		then
			sed '/^#/d' depth.csv > depth_clean.csv

		fi

		

		if [ ! -f hits.ls ]
		then
			awk '{print $3}' aligned.sam | \
			 sort | \
			 uniq -c | \
			 sort -nr >hits.ls 
			rm aligned.sam
		fi
	done
}

export -f bowtie_align
echo $READSPATH



parallel --jobs 15 bowtie_align ::: $READSPATH/*
