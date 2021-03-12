library(ggplot2)
library(dplyr)
library(ggsci)
library(gridExtra)
library(randomForest)
library(data.table)

# Load plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')
source("r_functions_and_parameters/r_options_BipEx.r")

GNOMAD_CHECK <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/19_gnomAD_check.tsv'
save_figures <- TRUE

df <- fread(GNOMAD_CHECK, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)

aes <- aes_string(x='not_inGnomAD_count', y='singleton_count', color='phenotype.PHENOTYPE_COARSE')
p <- create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'19_gnomAD_check_ColPheno'),
  n_x_ticks=5, x_label='Count of Singletons to in gnomAD', y_label='Count of singletons')
print(p)

aes <- aes_string(x='not_inGnomAD_count', y='singleton_count', color='phenotype.PROJECT_OR_COHORT')
p <- create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'19_gnomAD_check_ColBatch'),
  n_x_ticks=5, x_label='Count of Singletons to in gnomAD', y_label='Count of singletons')
print(p)

aes <- aes_string(x='not_inGnomAD_count', y='singleton_count', color='phenotype.LOCATION')
p <- create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'19_gnomAD_check_ColLocation'),
  n_x_ticks=5, x_label='Count of Singletons to in gnomAD', y_label='Count of singletons')
print(p)

df_SWE <- df %>% filter(phenotype.LOCATION == "Stockholm, SWE")
aes <- aes_string(x='not_inGnomAD_count', y='singleton_count', color='phenotype.PHENOTYPE_COARSE')
p <- create_pretty_scatter(df_SWE, aes, save_figure=FALSE, file=paste0(PLOTS,'19_gnomAD_check_ColPheno'),
  n_x_ticks=5, x_label='Count of Singletons to in gnomAD', y_label='Count of singletons')
print(p)