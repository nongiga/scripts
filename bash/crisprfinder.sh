
#!/bin/bash
# perl ${CASDIR}/CRISPRCasFinder.pl \
# -cas \
# -in assembly/Sample_Maccabi_Ecoli_SeqPlate9_H9/contigs.fasta \
# -gff prokka/Sample_Maccabi_Ecoli_SeqPlate9_H9/Sample_Maccabi_Ecoli_SeqPlate9_H9.gff \
# -faa prokka/Sample_Maccabi_Ecoli_SeqPlate9_H9/Sample_Maccabi_Ecoli_SeqPlate9_H9.faa \
# -out crisprfinder/Sample_Maccabi_Ecoli_SeqPlate9_H9 \
# -so ${CASDIR}/sel392v2.so



crisprfinder() {

	BASENAME=$(basename $1)

	
	if [ ! -f crisprfinder/$BASENAME/result.json ]
	then
		rm -rf ./crisprfinder/$BASENAME/result.json

		perl ${CASDIR}/CRISPRCasFinder.pl \
		-cas -rcfowce \
		-in $1/contigs.fasta \
		-gff prokka/$BASENAME/$BASENAME.gff \
		-faa prokka/$BASENAME/$BASENAME.faa \
		-out crisprfinder/$BASENAME \
		-so ${CASDIR}/sel392v2.so
	fi
}

export CASDIR="/media/kishonylab/KishonyStorage/Apps/CRISPRCasFINDER/CRISPRCasFinder-release-4.2.20"
export PATH=/media/kishonylab/KishonyStorage/Apps/CRISPRCasFINDER/CRISPRCasFinder-release-4.2.20/bin/:$PATH

export -f crisprfinder
echo $READSPATH
parallel --jobs 1 crisprfinder ::: assembly/* 