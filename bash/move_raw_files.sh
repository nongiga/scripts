#!/bin/bash
############################################################################
# ARRANGE FILES
# Olga Snitser, 20 May 2020
# Modified by Noga Aharony, June 7, 2020
############################################################################
# Paths to input/output files and folders


############################################################################
# Arrange files

FOLDER_KEY="Sample_Maccabi_Ecoli_SeqPlate"
RAW_FOLDERS=("/media/kishonylab/KishonyStorage/Illumina_raw/*/$FOLDER_KEY*"  "/media/kishonylab/KishonyStorage/Illumina_raw/*/Raw_data/$FOLDER_KEY*")
SAVE_PATH="/media/kishonylab/KishonyStorage/noga/MaccabiUTI/Filtered_data/"

mkdir -p $SAVE_PATH

move_raw_files() {

	EXT="_combined.fastq"

	FOLDER_NAME=$1
	SAVE_PATH=$2


	echo "copying over & combining reads from $FOLDER_NAME"
	SAVE_FOLDER_NAME="$SAVE_PATH$(basename $FOLDER_NAME)"
	mkdir -p $SAVE_FOLDER_NAME
	[ ! -f $SAVE_FOLDER_NAME/R1$EXT ] && zcat $FOLDER_NAME/*R1_001.fastq.gz > $SAVE_FOLDER_NAME/R1$EXT
	[ ! -f $SAVE_FOLDER_NAME/R2$EXT ] && zcat $FOLDER_NAME/*R2_001.fastq.gz > $SAVE_FOLDER_NAME/R2$EXT

}
export -f move_raw_files

 parallel --jobs 60 move_raw_files ::: ${RAW_FOLDERS[@]} ::: $SAVE_PATH