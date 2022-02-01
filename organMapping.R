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
  make_option(c('-a', '--annot_file'), type = 'character'
              ,default = NULL
              ,help = 'Annotation file name to collapse tissues to organs'
              ,metavar = 'character')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

binary_matrix <- read.delim(opt$file, header = T, check.names = F)
metadata <- read.delim(opt$annot_file, header = T, stringsAsFactors = T)

organCollapser <- function(consensus_elements, annotation_file) {

	niveles <- levels(annotation_file$Facet.ontology.ID)

	rv <- data.frame(chr = consensus_elements$chr
					 , start = consensus_elements$start
					 , end = consensus_elements$end)

	for (nivel in niveles) {
	  
	  annotation_file %>% filter(Facet.ontology.ID == nivel) -> tmp
	  tmp[,3] %>% as.character() %>% unique() %>% print() -> tmp
	  tmp_submatrix <- consensus_elements %>% select(all_of(tmp))
	  tmp_submatrix <- tmp_submatrix %>% rowSums(.)
	  tmp_submatrix <- ifelse(test = tmp_submatrix > 0, yes = 1, no = 0)
	  rv <- cbind(rv, tmp_submatrix)
	}
	
	names(rv) <- c("chr", "start", "end", niveles)

	return(rv)
}

collapsed_file <- organCollapser(binary_matrix, metadata)
output_name <- 'consensus_elements_annot_collapsed.txt'

write.table(collapsed_file, output_name
			, col.names =T, sep="\t"
			, row.names = F, quote=F)
