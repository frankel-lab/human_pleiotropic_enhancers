#! /bin/bash

#Author:laiker96

ENCODE_DNASE_METADATA=$1
ENCODE_H3K27ac_METADATA=$2
THREADS=$3
OVERLAP_FRACTION_REPLICATES=$4
OVERLAP_MERGE_FRACTION=$5
ENCODE_ANNOT_FILE=$6

#Find DNaseI-seq enriched regions with MACS2
call_dnase_peaks.sh $ENCODE_DNASE_METADATA $THREADS
#Intersect DNaseI-seq data with H3K27ac ChIP-seq data
intersect_dnase_h3k27ac.sh $ENCODE_H3K27ac_METADATA

#Create summits confidence intervals (gappedPeaks) for each tissue
for FILE in *_SNRs.txt; do
	N_LINES=$(wc -l < $FILE)
	echo $N_LINES
	if [[ $N_LINES -gt 1 ]]; then
		#If experiment has replicates create CI for summits
		make_dataframes.sh "$FILE" "$OVERLAP_FRACTION_REPLICATES"
		parallel_peaks.R -f "$FILE" -c "$THREADS"
	else
		FILE_ID=$(cut -f 1 $FILE)
		narrow2gapped.sh "$FILE_ID"_peaks.narrowPeak "$FILE_ID"
	fi
done

#Sort gappedPeak files
for FILE in *.gappedPeak; do

	bedtools sort -i "$FILE" \
		| awk 'BEGIN {FS="\t"; OFS="\t"} {if ($3 - $2 == 0) $3=$3+1; print$0}' \
		> tmp && mv tmp "$FILE"

done
#Create summit clusters file
cat *.gappedPeak | cut -f 1,2,3,4 | bedtools sort \
	| bedtools merge -c 4 -o distinct | sed 's/;/,/g' > summit_clusters.bed     

#Keep only autosomes and chrX
grep -wP "chr[0-9]+|chrX" summit_clusters.bed > tmp && mv tmp summit_clusters.bed
#Create consensus elements from summit clusters using MACS2 previously called peaks
create_consensus_elements.py summit_clusters.bed
#Resolve overlapping peaks
resolve_peaks.sh consensus_elements.bed "$OVERLAP_MERGE_FRACTION"
#Annotate peak activity

cp consensus_elements_resolved.bed consensus_elements_resolved_annot.bed

for FILE in *_SNRs.txt; do

	bedmap --echo --delim '\t'  --fraction-map 1.00 --indicator \
		consensus_elements_resolved_annot.bed H3K27ac"$(tail -n 1 "$FILE" | cut -f 1)"_peaks.narrowPeak \
		> tmp && mv tmp consensus_elements_resolved_annot.bed
		
	echo "$FILE" | awk -F "_" '{print $1}' >> header.txt
	
done

tr '\n' '\t' < header.txt | paste <(echo -e "chr\tstart\tend") - \
	| cat - consensus_elements_resolved_annot.bed > consensus_elements_resolved_annot.txt
#Cleanup
rm header.txt
#Collapse tissues to UBERON organs
organMapping.R -f consensus_elements_resolved_annot.txt -a "$ENCODE_ANNOT_FILE"
