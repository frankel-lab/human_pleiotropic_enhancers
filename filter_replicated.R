#! /usr/bin/env Rscript
#Author:laiker96
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
# Argument parser
option_list = list(
  make_option(c('-f', '--file'), type = 'character'
              ,default = NULL
              ,help = 'Input file name'
              ,metavar = 'character'),
  make_option(c('-o', '--outname'), type = 'character'
              ,default = NULL
              ,help = 'Output file name'
              ,metavar = 'character')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#Read H3K27ac peaks file
df <- read.table(file = opt$file, header = F, sep = "\t")
outname <- opt$outname


filterReplicatedPeaks <- function(df) {
	
	if (dim(df)[2] >= 5) {
	  
		active <- rowSums(df[, -c(1:3)])
		total <- (dim(df)[2] - 3)
		fraction <- (active/total)
		
		if (dim(df)[2] > 5) {
			df[fraction >= 0.5, c(1,2,3)] -> df
		}
		else {
			df[fraction > 0.5, c(1,2,3)] -> df
		}
	}
	return(df)
}

df <- filterReplicatedPeaks(df)
write.table(x = df, file = outname, quote = F, sep = "\t", row.names = F, col.names = F)
