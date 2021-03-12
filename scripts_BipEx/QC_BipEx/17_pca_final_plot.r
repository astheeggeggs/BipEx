library(ggplot2)
library(dplyr)
library(ggsci)
library(data.table)

source('r_functions_and_parameters/pretty_plotting.r')
source('r_functions_and_parameters/r_options_BipEx.r')

FINAL_SAMPLE_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz | gzcat'

df <- fread(FINAL_SAMPLE_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df)),]
save_figures <- TRUE

for (i in c(1,3,5)) {
  aes <- aes_string(x=paste0('PCA.PC',i), y=paste0('PCA.PC',i+1), color='PHENOTYPE_COARSE')
  create_pretty_scatter(df, aes,
  	save_figure=save_figures,
  	file=paste0(PLOTS,'17_PC',i,'_PC',i+1,'_final_PCs'),
  	x_label=paste0('Principal Component ', i),
  	y_label=paste0('Principal Component ', i+1))
  aes <- aes_string(x=paste0('PCA.PC',i), y=paste0('PCA.PC',i+1), color='LOCATION')
  create_pretty_scatter(df, aes,
  	save_figure=save_figures,
  	file=paste0(PLOTS,'17_PC',i,'_PC',i+1,'_final_PCs_Location'),
  	x_label=paste0('Principal Component ', i),
  	y_label=paste0('Principal Component ', i+1))
}
