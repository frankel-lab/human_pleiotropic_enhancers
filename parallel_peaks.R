#! /usr/bin/env Rscript
#Author:laiker96
suppressMessages(library(doParallel))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
# Argument parser
option_list = list(
  make_option(c('-f', '--file'), type = 'character'
              ,default = NULL
              ,help = 'Text file with file ids and SNRs'
              ,metavar = 'character'),
  make_option(c('-c', '--cores'), type = 'integer'
              ,default = "1"
              ,help = 'Number of cores for parallel computing'
              ,metavar = 'character')
  
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

#######################################
list_dataframes <- read.table(opt$file, sep = '\t', header = F, dec = ',')
higher_SNR_file <- paste('H3K27ac', tail(list_dataframes, 1)[1], sep='')
list_dataframes <- list.files(path = ".", pattern = paste(higher_SNR_file, "_ENCF.*.bed", sep=''))
print(list_dataframes)



list_to_join <- list()

for (j in list_dataframes) {
  
  tmp <- read.table(j, header = F, sep = '\t', stringsAsFactors = TRUE)
  assign(x = strsplit(x = j, split = "\\.")[[1]][1], value = tmp)
  list_to_join[[j]] <- tmp
  rm(tmp)
  
}

df <- do.call(what = rbind, args = list_to_join)
df %>% arrange(V4) -> df

make_ci <- function(peaks) {
  clusters_n <- length(levels(peaks[, 4])) # Set output dataframe row length
  # Names according to the gappedPeak format
  col_names <- c("chrom", "chromStart", "chromEnd", "name", "score"
                 , "strand", "thickStart", "thickEnd"
                 , "itemRgb", "blockCount", "blockSizes"
                 , "blockStarts", "signalValue", "pValue"
                 , "qValue")
  #Create empty output dataframe
  rv <- data.frame(matrix(nrow = clusters_n, ncol = 15, dimnames = list(NULL, col_names)))
  peak_elements <- levels(peaks[, 4])
  j <- 1
  for (peak_element in peak_elements) {
    df_to_append <- data.frame(matrix(nrow = 1, ncol = 15, dimnames = list(NULL, col_names)))
    peaks %>% filter(V4 == peak_element) -> tmp
    #df_to_append["name"] <- c(as.character(tmp[1, 4]), as.character(tmp[, 14])) %>% paste(., collapse = ":")
    cluster_name <- c(paste(tmp[1,1], tmp[1,2], tmp[1,3], sep = "-"), paste(tmp[,11], tmp[,12], tmp[,13], sep = "-"))
    cluster_name <- paste(cluster_name, collapse = ";")
    df_to_append["name"] <- cluster_name
    df_to_append["chrom"] <- as.character(tmp[1, 1])
    df_to_append["thickStart"] <- 0
    df_to_append["thickEnd"] <- 0
    df_to_append["itemRgb"] <- "0,0,0"
    cluster_boundaries <- c(tmp[1, 2] + tmp[1, 10], tmp[, 12] + tmp[, 20])
    minima <- min(cluster_boundaries)
    maxima <- max(cluster_boundaries)
    if (minima == maxima) {
      maxima <- maxima + 1
    }
    df_to_append["chromStart"] <- minima
    df_to_append["chromEnd"] <- maxima
    df_to_append["score"] <- mean(c(tmp[1, 5], tmp[, 15]))
    df_to_append["strand"] <- "."
    df_to_append["blockCount"] <- length(unique(cluster_boundaries))
    df_to_append["blockSizes"] <- paste(rep(1, df_to_append["blockCount"]), collapse = ",")
    df_to_append["blockStarts"] <- paste(unique(cluster_boundaries - as.integer(df_to_append["chromStart"])), collapse = ",")
    df_to_append["signalValue"] <- mean(c(tmp[1, 7], tmp[, 17]))
    # Combine P values according with Fisher's Method
    pValues <- c(tmp[1, 8], tmp[, 18])
    fisher_X <- 2/log(10) * sum(pValues)
    fisher_p <- pchisq(fisher_X, df = length(pValues) * 2, lower.tail = F)
    df_to_append["pValue"] <- -log(fisher_p)
    # Q value undefined
    df_to_append["qValue"] <- (-1)
    rv[j, ] <- df_to_append
    j <- j + 1

  }
  
  return(rv)
  
  }


chunk2 <- function(x,n) {
	split(x, cut(seq_along(x), n, labels = FALSE)) 
}


split_df <- function(df, clusters) {

  #Split dataframe into chunks for parallel forking processing
  rv <- list()
  peak_elements <- levels(df[, 4])
  chunks <- chunk2(x = peak_elements, n = clusters) 
  i <- 1
  
  for (chunk in chunks) {
    df %>% filter(V4 %in% chunk) -> tmp
    tmp[, 4] <- factor(tmp[, 4])
    rv[[i]] <- tmp
    i <- i + 1
  }
  return(rv)
}

dfs <- split_df(df = df, clusters = opt$cores)
cl <- parallel::makeForkCluster(opt$cores)
doParallel::registerDoParallel(cl)
gp <-foreach(i=1:opt$cores, .combine=rbind) %dopar%
  make_ci(peaks = dfs[[i]])

output_name <- strsplit(list_dataframes[1], "_")[[1]][1]
file_name = paste(output_name, "gappedPeak", sep = ".")

write.table(x = gp, file = file_name, quote = F
			, row.names = F, col.names = F, sep = "\t")
