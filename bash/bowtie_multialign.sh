#!/bin/bash
############################################################################
# multialign to bowtie + convert to # of hits
# Noga Aharony
# 01 July 2020
############################################################################


bowtie_multialign() {
	read -r subdir bp v <<<$(echo "$1 $2 $3")
	for SP in $subdir/*.gff ; do
		SEQDIR=$(basename -- $SP)
		SEQDIR="${SEQDIR%.*}"

		
		cd $subdir/$SEQDIR
		if [ ! -f aligned_sorted$bp.multi.bam ]
		then 
			bowtie \
			-S -p 1 -v $v -a \
			--un unaligned$bp.multi.fastq \
			$1/pangenome \
			R1_spl$bp.combined.trimmed.fastq.gz \
			aligned$bp.multi.sam #outputq

			samtools view -bS -o aligned$bp.multi.bam aligned$bp.multi.sam
			samtools sort  aligned$bp.multi.bam -o aligned_sorted$bp.multi.bam
			samtools index aligned_sorted$bp.multi.bam
			rm aligned$bp.multi.bam aligned$bp.multi.sam
		fi
		
		if [ ! -f aligned$bp.multi.vcf ]
		then
			samtools mpileup -t AD -ugf $subdir/roary_output/pan_genome_reference.fa aligned_sorted$bp.multi.bam > pileup$bp.multi.vcf.temp
			/media/kishonylab/KishonyStorage/Apps/bcftools/bcftools/bcftools view -o aligned$bp.multi.vcf pileup$bp.multi.vcf.temp
			rm pileup$bp.multi.vcf.temp
		fi

		if [ ! -f depth_clean$bp.multi.csv ]
		then
			awk -F'[=\t;]' '{print $1 , $2, $9}' aligned$bp.multi.vcf | sed '/^#/d' > depth_clean$bp.multi.csv
		fi

	done
}


export -f bowtie_multialign


read -r maindir bp threads <<<$(echo "$1 $2 $3")
v=1
if [ $bp -gt 40 ]
then
	v=3
fi

parallel --jobs $threads bowtie_multialign ::: $maindir/* ::: $bp ::: $v
