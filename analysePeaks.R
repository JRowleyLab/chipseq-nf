# Store arguments of:
# file
# fold change
args <- commandArgs(trailingOnly = TRUE)

xls <- args[1]
fold_change <- as.numeric(args[2])

#libraries
library(tidyverse)

# Read in .xls peaks file
diffpeaks <- read.csv(xls, sep = "\t")


diffpeaks %>% 
  mutate(log2fc = log2(.[[9]] / .[[10]]) ) %>%
  filter(log2fc < -fold_change | log2fc > fold_change) -> diffpeaks_out


# Write filtered peaks files
name = paste0("diffpeaks_",fold_change,"fc.tsv")

write.table(diffpeaks_out, file = name, sep = "\t", quote = F, row.names = F)
