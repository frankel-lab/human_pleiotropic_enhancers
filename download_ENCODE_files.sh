#! /bin/bash

#Author:laiker96
#Dependencies:{gnu-parallel, curl, samtools}
#The ENCODE_METADATA_FILE is a tab delimited file downloaded from the 
#ENCODE Portal with the file accession name, experiment accession,
#tissue or cell type and other metadata. See H3K27ac_datasets_ENCODE_Roadmaps_metadata.tsv
#or DNase_datasets_ENCODE_Roadmaps_metadata.tsv for examples. It also contains the
#ENCODE Portal links to all the files.
ENCODE_METADATA_FILE=$1 #tsv file with ENCODE metadata as described before
TYPE_OF_FILE=$2 #BAM of narrowPeak
THREADS=$3 #Number of threads to use for parallel downloading

function DownloadEncodeData {
	
	
	ID=$(basename "$1" .tsv)
	FILE_DIR=$(dirname "$1")
	
	#Extract file URLs
	cut -f 47 "$1" | tail -n +2 > "$FILE_DIR"/"$ID".links

	if [[ "$2" == "BAM" ]]; then
	#Files are alignments in .bam format

		parallel -a "$ID."links -j "$3" curl -O -J -L {}
		
		#Index the BAM files
		for BAM in *.bam; do
			samtools index -@ "$3" "$BAM"
		done

	#Files in narrowPeak (bed.gz) format
	else

		parallel -a "$ID."links -j "$3" curl -O -J -L {}
	fi

	#Cleanup files
	rm "$ID."links

}

DownloadEncodeData $ENCODE_METADATA_FILE $TYPE_OF_FILE $THREADS
