library(ggplot2)
library(dplyr)
# library(ggsci)
library(data.table)

source('r_functions_and_parameters/pretty_plotting.r')
source("r_functions_and_parameters/r_options_BipEx.r")
source("r_functions_and_parameters/helpful_functions.r")

PCA_EUR_SCORES <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/11_pca_scores.strict_european.tsv"
URV_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv'

save_figures <- TRUE

df_PCs <- fread(PCA_EUR_SCORES, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df_URVs <- fread(URV_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df_URVs)[names(df_URVs) == 'Samples'] <- 'sample'
df <- merge(df_PCs, df_URVs)

for (i in c(1,3,5)) {
    aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.LOCATION')
    p <- create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'11_PC',i,'_PC',i+1,'_EUR'),
        x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    # p <- p + geom_point(data=df %>% filter(keep == FALSE),
    #   mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
    #   inherit.aes=FALSE, shape=1, show.legend=FALSE)
    print(p)
}
