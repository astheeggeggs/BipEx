library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

# Get thresholds and file locations.
source("r_functions_and_parameters/r_options_BipEx.r")

IBD_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.tsv"
df <- fread(IBD_FILE, sep='\t', stringsAsFactors=TRUE, header=TRUE, data.table=FALSE)

df$inferred_relationship <- 'Siblings'
df$inferred_relationship[df$ibd.PI_HAT <= IBD_THRESHOLD] <- 'Unrelated'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 < 0.05] <- 'Duplicate/Monozygotic twins'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 > 0.9] <- 'Parent-Offspring'

df_related <- filter(df, inferred_relationship != 'Unrelated') %>% 
  mutate(i=droplevels(i), j=droplevels(j))

samples <- names(sort(table(unlist(df_related[,c('i', 'j')])), decreasing=TRUE))
removed <- character()

for (row in 1:nrow(df_related)) {
  
  i_sample <- as.character(df_related[row, 'i'])
  j_sample <- as.character(df_related[row, 'j'])
  
  i_index <- match(i_sample, samples)
  j_index <- match(j_sample, samples)
  
  if (is.na(i_index) | is.na(j_index)) { next }
  
  if (i_index <= j_index) {
    removed <- c(removed, i_sample)
    samples <- samples[samples != i_sample]
  } else {
    removed <- c(removed, j_sample)
    samples <- samples[samples != j_sample]
  } 
  
}

fwrite(as.data.frame(removed), file=SAMPLE_LIST_IBD, quote=FALSE, row.names=FALSE, col.names=FALSE)
