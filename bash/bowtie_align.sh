#!/bin/bash
############################################################################
# align to bowtie + convert to final 
# Noga Aharony
# 01 July 2020
############################################################################


bowtie_align() {
	read -r subdir bp v <<<$(echo "$1 $2 $3")
	for SP in $subdir/*.gff ; do
		SEQDIR=$(basename -- $SP)
		SEQDIR="${SEQDIR%.*}"

		
		cd $subdir/$SEQDIR

		if [ ! -f aligned_sorted$bp.bam ]
		then 
			bowtie \
			-S -p 1 -v $v -m 1 \
			--max multialigned$bp.fastq \
			--un unaligned$bp.fastq \
			$1/pangenome \
			R1_spl$bp.combined.trimmed.fastq.gz \
			aligned$bp.sam #outputq

			samtools view -bS -o aligned$bp.bam aligned$bp.sam
			samtools sort  aligned$bp.bam -o aligned_sorted$bp.bam
			samtools index aligned_sorted$bp.bam

			gzip multialigned$bp.fastq unaligned$bp.fastq

			rm aligned$bp.bam aligned$bp.sam R1_spl$bp.combined.trimmed.fastq.gz
		fi
		
		

		if [ ! -f depth_clean$bp.csv ]
		then
			samtools mpileup -t AD -ugf $subdir/roary_output/pan_genome_reference.fa aligned_sorted$bp.bam > pileup$bp.vcf.temp
			/media/kishonylab/KishonyStorage/Apps/bcftools/bcftools/bcftools view -o aligned$bp.vcf pileup$bp.vcf.temp
			awk -F'[=\t;]' '{print $1 , $2, $9}' aligned$bp.vcf | sed '/^#/d' > depth_clean$bp.csv
			rm pileup$bp.vcf.temp aligned$bp.vcf

		fi

	done
}


export -f bowtie_align


read -r maindir bp threads <<<$(echo "$1 $2 $3")
v=1
if [ $bp -gt 40 ]
then
	v=3
fi
echo "bowtie_align"
parallel --jobs $threads bowtie_align ::: $maindir/* ::: $bp ::: $v
