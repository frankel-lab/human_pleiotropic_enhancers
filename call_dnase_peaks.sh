#! /bin/bash

#Author:laiker96
#Dependencies: {macs2, gnu-parallel, bedops, samtools}
#ENCODE_METADATA_FILE is a tsv downloaded from ENCODE portal
#that contains information of each experiment, such as 
#accession numbers and UBERON ontology terms.

ENCODE_METADATA_FILE=$1
THREADS=$2
UBERON_IDS=$(cut -f 9 "$ENCODE_METADATA_FILE" | tail -n +2 | sort -g | uniq)


function checkExperimentParity {
	# $1:$ENCODE_METADATA_FILE; $2:$THREADS
	while read k; do
		#Extract ENCODE accession numebers
		BAM_FILES=$(grep -wf <(echo "$k") "$1" | cut -f 1)
		
		while read bam; do

			FILE=$(echo "$bam".bam)
			INDEX=$(echo "$bam".bam.bai)
			#Check parity
			N=$(samtools view --threads "$2" -c -F 1 "$FILE") 
			
			if [[ $N -eq 0  ]]; then

				echo "$FILE" >> PE_bamFiles.txt
				echo "$k" >> PE_UBERON_ids.txt
			else
				echo "$FILE" >> SE_bamFiles.txt
				echo "$k" >> SE_UBERON_ids.txt
			fi


		done < <(echo "$BAM_FILES")


	done < <(echo "$UBERON_IDS")

}

function computeSNR {

	FIP=$(bedmap --echo --delim '\t' --fraction-map 1 --count \
		"$1" <(bam2bed --max-mem 8G < "$2") \
		| awk '{i+=$11} END{print i}')

	
	TF=$(bam2bed --max-mem 8G < "$2" | wc -l)

	echo "scale=3; $FIP/$TF" | bc
}

function callParallelPeaksPE {
	
	#PE_bamFiles:$1; THREADS:$2; FILE:$3
	parallel -a "$1" -j "$2" bash call_peaks.sh {} {//} {/.} "PE"

	while read LINE; do

		FILE=$(echo "$LINE" | cut -f 1)
		UBERON_ID=$(echo "$LINE" | cut -f 2)
		ID=$(echo "$FILE" | awk -F "." '{print $1}')
		BASE_NAME=$(basename "$FILE" .bam)
	
		SNR=$(computeSNR "$ID"_peaks.narrowPeak "$FILE")
		paste <(echo $BASE_NAME) <(echo $SNR) >> "$UBERON_ID"_SNRs.txt
		sed 's/\./,/g' "$UBERON_ID"_SNRs.txt | sort -gk2 > tmp
		mv tmp "$UBERON_ID"_SNRs.txt

	done < <(paste "$1" "$3")
}	


function callParallelPeaksSE {
	
	#SE_bamFiles:$1; THREADS:$2; FILE:$3
	parallel -a "$1" -j "$2" bash call_peaks.sh {} {//} {/.} "SE"

	while read LINE; do
		
		FILE=$(echo "$LINE" | cut -f 1)
		UBERON_ID=$(echo "$LINE" | cut -f 2)
		ID=$(echo "$FILE" | awk -F "." '{print $1}')
		BASE_NAME=$(basename "$FILE" .bam)
	
		SNR=$(computeSNR "$ID"_peaks.narrowPeak "$FILE")
		paste <(echo $BASE_NAME) <(echo $SNR) >> "$UBERON_ID"_SNRs.txt
		sed 's/\./,/g' "$UBERON_ID"_SNRs.txt | sort -gk2 > tmp
		mv tmp "$UBERON_ID"_SNRs.txt
	
	done < <(paste "$1" "$3")
	
}

checkExperimentParity "$ENCODE_METADATA_FILE" "$THREADS"

if [ -f PE_bamFiles.txt ]; then
	callParallelPeaksPE PE_bamFiles.txt "$THREADS" PE_UBERON_ids.txt
fi

if [ -f SE_bamFiles.txt ]; then
	callParallelPeaksSE SE_bamFiles.txt "$THREADS" SE_UBERON_ids.txt
fi
