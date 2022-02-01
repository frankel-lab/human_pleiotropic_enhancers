#! /bin/bash

#Author:laiker96
#Convert ENCODE narrowPeak format to ENCODE gappedPeak format
#Format specifications = {

#gappedPeak = https://genome.ucsc.edu/FAQ/FAQformat.html#format14
#narrowPeak = https://genome.ucsc.edu/FAQ/FAQformat.html#format12

#}

PEAK=$1
PEAK_ID=$2
function narrowToGapped {

	awk '{print $1 "\t" $2 + $10 "\t" $2 + $10 + 1 "\t" $1"-"$2"-"$3}' "$1" \
		> H3K27ac"$2".gappedPeak

}

narrowToGapped "$PEAK" "$PEAK_ID"
