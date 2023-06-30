# Store arguments of:
# file
# fold change
args <- commandArgs(trailingOnly = TRUE)

xls <- args[1]
fold_change <- as.numeric(args[2])
outdir <- args[3]


#libraries
library(tidyverse)

# Read in .xls peaks file
diffpeaks <- read.csv(xls, sep = "\t")


diffpeaks %>% 
  mutate(log2fc = log2(.[[9]] / .[[10]]) ) %>%
  filter(log2fc < -fold_change | log2fc > fold_change) -> diffpeaks_out


# Write filtered peaks files
sample_name = gsub(".xls", "", basename(xls))
name = paste0(sample_name, "_diffpeaks_",fold_change,"log2fc.tsv")

write.table(diffpeaks_out, file = paste0(outdir, "/", name), sep = "\t", quote = F, row.names = F) # nolint
