#! /bin/bash

#Author:laiker96
#Script for peak calling DNaseI-seq data. 
#Input:BAM file
#Output:narrowPeak file

FILE=$1 #BAM file name
OUTDIR=$2 #Output directory name
OUTNAME=$3 #Output name prefix
PARITY=$4 #Parity(PE or SE)


function CallPeaksSE {
	#Call peaks for SE experiments	
	macs2 callpeak -t "$1" -f BAM -g hs --outdir "$2" -n "$3" \
		--nomodel --extsize 200 -q 0.01 --min-length 50

}

function CallPeaksPE {
	#Call peaks for PE experiments	
	macs2 callpeak -t "$1" -f BAMPE -g hs --outdir "$2" -n "$3" \
		-q 0.01 --min-length 50

}


#Call functions
if [[ "$PARITY" == "SE" ]]; then

	CallPeaksSE "$FILE" "$OUTDIR" "$OUTNAME"
else

	CallPeaksPE "$FILE" "$OUTDIR" "$OUTNAME"

fi
