#! /bin/bash

#Author:laiker96
#Dependencies: {bedtools, R}
#ENCODE_METADATA_FILE is a tsv downloaded from ENCODE portal
#that contains information of each experiment, such as 
#accession numbers and UBERON ontology terms

ENCODE_METADATA_FILE=$1


function SortBedFiles {

	for BED in *.bed; do
		bedtools sort -i "$BED" > tmp
		mv tmp "$BED"
	done

}

function intersectH3K27acDNaseData {

	UBERON_IDS=$(cut -f 9 "$1" | tail -n +2 | sort -g | uniq)

	while read k; do

		BED_FILES=$(grep -w $k "$1" | cut -f 1)
		
		while read bed; do

			FILE=$(echo "$bed".bed)
			cat "$FILE" >> "$k"_mergedPeaks.bed

		done < <(echo "$BED_FILES")
			
		bedtools sort -i "$k"_mergedPeaks.bed | bedtools merge > tmp
		mv tmp "$k"_mergedPeaks.bed

		while read bed; do

			FILE=$(echo "$bed".bed)
			bedmap --echo --delim "\t" --indicator "$k"_mergedPeaks.bed \
				"$FILE" > tmp
			mv tmp "$k"_mergedPeaks.bed

		done < <(echo "$BED_FILES")
		
		filter_replicated.R -f "$k"_mergedPeaks.bed -o filtered_"$k"_mergedPeaks.bed
		
		bedtools intersect -wa -a "$(cut -f 1 "$k"_SNRs.txt | tail -n 1)"_peaks.narrowPeak \
			-b filtered_"$k"_mergedPeaks.bed \
			| bedtools sort \
			| uniq > H3K27ac"$(cut -f 1 "$k"_SNRs.txt | tail -n 1)"_peaks.narrowPeak
		

	done < <(echo "$UBERON_IDS")
}

#Sort all bed files in the current directory
SortBedFiles 
#Process peak data
intersectH3K27acDNaseData $ENCODE_METADATA_FILE
