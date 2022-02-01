#! /bin/bash

#Author:laiker96
#Arguments description:
#SNR_FILE: String
#Narrowpeak file (10 fields) of the highest signal-to-noise 
#ratio DNase-seq Alignment
SNR_FILE=$1
#FRACTION_OVERLAP: Float 
#Fraction of overlap between peaks in the HIGHER_SNR_FILE and the other
#replicates files. Not reciprocal (the fraction of overlap must be fulfilled
#by either the HIGHER_SNR_FILE or the other replicates)
FRACTION_OVERLAP=$2

#This function takes two narrowPeaks ($1 and $2) and one float ($3) as 
#arguments and writes a 20-field dataframe of the 10 fields of each
#peak file that overlap each other by the fraction or more specified in $3.
#If there are mutiple overlaps, it keeps the one in which the summit position
#is smaller between the peaks.
function make_dataframes {

	#Create Dataframe with bedmap
	bedmap --echo --skip-unmapped --fraction-either $3 --echo-map "$1" "$2" > tmp
	#Split the dataframe into multiple overlap peaks and single overlap peaks
	grep ";" tmp > tmp_multiple
	grep -v ";" tmp > tmp_single

	#Keep the overlap based on the lower distance between summits
	while read p; do

		PEAK=$(echo "$p" | awk -F "|" '{print $1}')
		REPLICATES=$(echo "$p" | awk -F "|" '{print $2}')
		REPLICATES=$(tr ';' '\n' < <(echo "$REPLICATES"))
		SUMMIT_POSITION=$(awk '{print $2 + $10}' < <(echo "$PEAK"))
		REPLICATE_SUMMIT_POSITION=$(awk '{print $2 + $10}' < <(echo "$REPLICATES" | head -n 1))
		CURRENT_DELTA=$(expr $SUMMIT_POSITION - $REPLICATE_SUMMIT_POSITION | sed 's/-//')
		RV=$(echo "$REPLICATES" | head -n 1)

		while read j; do

			SUMMIT_REPLICATE=$(awk '{print $2 + $10}' < <(echo "$j"))
			DELTA=$(expr "$SUMMIT_POSITION" - "$SUMMIT_REPLICATE" | sed 's/-//')
			if [[ "$DELTA" -le "$CURRENT_DELTA" ]]; then
				CURRENT_SUMMIT=$DELTA
				RV=$j
			
			fi
				

		done < <(echo "$REPLICATES" | tail -n +2)
		echo "$RV" >> tmp_to_paste

	done < tmp_multiple

	#Create file prefixes
	PEAK_ID=$(echo "$1" | awk -F "_" '{print $1}')
	REPLICATE_ID=$(echo "$2" | awk -F "_" '{print $1}')
	#Create single overlap file
	awk -F "|" '{print $1}' tmp_multiple > tmp_multiple_coords

	#Concatenate the single overlap files into a single file
	cat tmp_single <(paste -d "|" tmp_multiple_coords tmp_to_paste) | bedtools sort \
		| tr '|' '\t' > "$PEAK_ID"_"$REPLICATE_ID".bed
	
	#Clean up
	rm tmp*

}


HIGHER_SNR_FILE=$(cut -f 1 "$SNR_FILE" | tail -n 1)
REPLICATES=$(grep -v "$HIGHER_SNR_FILE" "$SNR_FILE" | cut -f 1)
HIGHER_SNR_FILE="H3K27ac"$HIGHER_SNR_FILE"_peaks.narrowPeak"

while read LINE; do
	#narrowPeak file name
	i="$LINE"_peaks.narrowPeak
	echo "Processing files..."
	echo "$HIGHER_SNR_FILE" "<---->" "$i"
	make_dataframes $HIGHER_SNR_FILE "$i" $FRACTION_OVERLAP
	make_dataframes "$i" $HIGHER_SNR_FILE $FRACTION_OVERLAP


	PEAK_ID=$(echo "$HIGHER_SNR_FILE" | awk -F "_" '{print $1}')
	REPLICATE_ID=$(echo "$i" | awk -F "_" '{print $1}')

	paste <(cut -f 1-10 --complement "$REPLICATE_ID"_"$PEAK_ID".bed) <(cut -f 1-10 \
		"$REPLICATE_ID"_"$PEAK_ID".bed) > reverted_dataFrame.bed

	mv reverted_dataFrame.bed "$REPLICATE_ID"_"$PEAK_ID".bed

	parallel --pipe --block 10M --ungroup LC_ALL=C grep -wFf "$REPLICATE_ID"_"$PEAK_ID".bed  \
		< "$PEAK_ID"_"$REPLICATE_ID".bed > dataFrame.txt
	
	rm "$REPLICATE_ID"_"$PEAK_ID".bed
	mv dataFrame.txt "$PEAK_ID"_"$REPLICATE_ID".bed
	

done < <(echo "$REPLICATES")
