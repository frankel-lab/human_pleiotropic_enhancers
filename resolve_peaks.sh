#! /bin/bash

FILE=$1
T=$2

#First sort bed file
bedtools sort -i $FILE > tmp && mv tmp $FILE


function resolvePeaks {

	bedmap --echo --delim ":" --fraction-either $2 --count $1 \
		| awk -F ":" '$2 == 1 {print $1}' \
		| cut -f 1,2,3 | bedtools sort > non_overlapping_peaks.bed


	bedmap --echo --delim ":" --fraction-either $2 --count $1 \
		| awk -F ":" '$2 != 1 {print $1}' \
		| cut -f 1,2,3 | bedtools sort > overlapping_peaks.bed

	while [[ -s overlapping_peaks.bed ]]; do
		
		bedmap --echo-map-range --fraction-either $2 overlapping_peaks.bed \
			| bedtools sort | uniq  > tmp
		cat tmp non_overlapping_peaks.bed | bedtools sort | uniq > peaks.bed

		bedmap --echo --delim ":" --fraction-either $2 --count peaks.bed \
			| awk -F ":" '$2 == 1 {print $1}' > non_overlapping_peaks.bed
		
		bedmap --echo --delim ":" --fraction-either $2 --count peaks.bed \
			| awk -F ":" '$2 != 1 {print $1}' > overlapping_peaks.bed


		done

	#Cleanup
	rm tmp overlapping_peaks.bed peaks.bed
	mv non_overlapping_peaks.bed consensus_elements_resolved.bed

}

resolvePeaks $FILE $T 
