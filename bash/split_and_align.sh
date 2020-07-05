#!/bin/bash
############################################################################
# splits sequnece to $2bp and 
# Olga Snitser
# 04 June 2020
# Modified by Noga Aharony, June 8, 2020
############################################################################
# Paths to input/output files and folders
READSPATH=$1

############################################################################
# Run roary
# $1: main dir, $2: bp length, $3: parallel, $4=input sequence folder, $5=multialigned

split_sequence(){
		for SP in $1/*.gff ; do
			SEQDIR=$(basename -- $SP)
			SEQDIR="${SEQDIR%.*}"

			SEQPATH=$3/$SEQDIR/
			mkdir -p $1/$SEQDIR; cd $1/$SEQDIR
			if [ ! -f R1_spl$2.combined.trimmed.fastq.gz ] && [ -f $SEQPATH/R1_combined.trimmed.fastq.gz ]
			then
				~/bin/seqkit sliding \
				-s $2 -W $2 \
				$SEQPATH/R1_combined.trimmed.fastq.gz | \
				gzip -c >R1_spl$2.combined.trimmed.fastq.gz
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

		v=1
		if [ $2 > 40 ]
		then
			v=3
		fi

		echo "$SEQDIR"
		cd $1/$SEQDIR

		$m=''
		echo $4
		if [ ! -f aligned$2.sam ] && [ $3 == 0 ]
		then 
			bowtie \
			-S -p 1 -v $v -m 1 \
			--max multialigned$2.fastq \
			--un unaligned$2.fastq \
			$1/pangenome \
			R1_spl$2.combined.trimmed.fastq.gz \
			aligned$2.sam #outputq
		fi

		if [ ! -f aligned$2_m.sam ] && [ $3 == 1 ]
		then 
			$m='_m'
			bowtie \
			-S -p 1 -v $v -a \
			--un unaligned$2$m.fastq \
			$1/pangenome \
			R1_spl$2$m.combined.trimmed.fastq.gz \
			aligned$2$m.sam #outputq
			
		fi
		if [ ! -f aligned_sorted$2$m.bam ]
		then 
			samtools view -bS -o aligned$2$m.bam aligned$2$m.sam
			samtools sort  aligned$2$m.bam -o aligned_sorted$2$m.bam
			samtools index aligned_sorted$2$m.bam
			rm aligned$2$m.bam aligned$2$m.sam
		fi
		
		if [ ! -f aligned$2$m.vcf ]
		then
			samtools mpileup -t AD -ugf $1/roary_output/pan_genome_reference.fa aligned_sorted$2$m.bam > pileup$2$m.vcf.temp
			/media/kishonylab/KishonyStorage/Apps/bcftools/bcftools/bcftools view -o aligned$2$m.vcf pileup$2$m.vcf.temp
			# another line here
			rm pileup$2$m.vcf.temp
		fi

		if [ ! -f depth_clean$2$m.csv ]
		then
			awk -F'[=\t;]' '{print $1 , $2, $9}' aligned$2$m.vcf | sed '/^#/d' > depth_clean$2$m.csv

		fi
	done
}
export -f bowtie_align


echo "split sequence"
parallel --jobs $3 split_sequence ::: $1/* ::: $2 ::: $4
echo "bowtie_build"
parallel --jobs $3 bowtie_build ::: $1/*
echo "bowtie_align"
parallel --jobs $3 bowtie_align ::: $1/* ::: $2 :: $5