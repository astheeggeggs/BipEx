library(data.table)
# Plotting functions:
source('QC_BipEx/r_functions_and_parameters/pretty_plotting.r')

PADDING_METRICS_FILE <- "gsutil cat gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/padding_sample_metrics.tsv.gz | gzcat"
df <- fread(PADDING_METRICS_FILE, stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE)

# Let's rename the PI column so it's not stupid.
save_figures <- FALSE

# Remove the low coverage samples as defined by Excel spreadsheet e-mailed from Laura Gauthier
low_coverage <- fread("QC_BipEx/Dalio_Low_Coverage_Samples.txt")$SAMPLE

df <- df[-which(df$s %in% low_coverage),]
df <- df[sample(nrow(df), replace=FALSE),]
PLOTS <- "padding_plots/"

width <- 160 
height <- 90
scaling <- 1
# CDFs, no splitting.
p <- create_pretty_cumulative(df, aes(qc_100.call_rate), 'Call Rate', 0,
    xlim=c(0.5,0.8), title='Call Rate', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_50.call_rate), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.call_rate), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'callRate_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_cumulative(df, aes(qc_100.dp_stats.mean), 'Mean Depth', 0,
    xlim=c(20, 70), title='Mean Depth', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_50.dp_stats.mean), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.dp_stats.mean), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'dpMean_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_cumulative(df, aes(qc_100.gq_stats.mean), 'Mean Genotype Quality', 0,
    xlim=c(35, 50), title='Mean Genotype Quality', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_50.gq_stats.mean), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.gq_stats.mean), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'gqMean_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')

df_50 <- df[,c('s', names(df)[grep("_50", names(df))])]
df_50$padding <- '50'
names(df_50) <- gsub("_50", "", names(df_50))
df_100 <- df[,c('s', names(df)[grep("_100", names(df))])]
names(df_100) <- gsub("_100", "", names(df_100))
df_100$padding <- '100'
df_0 <- df[,c('s', names(df)[grep("_0", names(df))])]
names(df_0) <- gsub("_0", "", names(df_0))
df_0$padding <- '0'

df_full <- bind_rows(df_0, df_50, df_100)
df_full$padding <- ordered(df_full$padding, levels = c("0", "50", "100"))

save_figures <- TRUE

create_pretty_boxplots(df_full, aes(x=padding, y=qc.n_called), aes(color=padding),
    x_label='n called', y_label='padding (bp)', title="", save_figure=save_figures,
    file=paste0(PLOTS,'n_called'), n_ticks=5)

create_pretty_boxplots(df_full, aes(x=padding, y=qc.r_ti_tv), aes(color=padding),
    x_label='Ti/Tv', y_label='padding (bp)', title="", save_figure=save_figures,
    file=paste0(PLOTS,'ti_tv'), n_ticks=5)

PADDING_METRICS_FILE <- "gsutil cat gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/padding_sample_metrics_ICE.tsv.gz | gzcat"
df <- fread(PADDING_METRICS_FILE, stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE)

# Remove the low coverage samples as defined by Excel spreadsheet e-mailed from Laura Gauthier
low_coverage <- fread("QC_BipEx/Dalio_Low_Coverage_Samples.txt")$SAMPLE

df <- df[-which(df$s %in% low_coverage),]
df <- df[sample(nrow(df), replace=FALSE),]
PLOTS <- "padding_plots/"

width <- 160 
height <- 90
scaling <- 1
# CDFs, no splitting.
p <- create_pretty_cumulative(df, aes(qc_150.call_rate), 'Call Rate', 0,
    xlim=c(0.9,1), title='Call Rate', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_100.call_rate), geom='line', pad=FALSE, col='grey30')
p <- p + stat_ecdf(aes(qc_50.call_rate), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.call_rate), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'callRate_ICE_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_cumulative(df, aes(qc_150.dp_stats.mean), 'Mean Depth', 0,
    xlim=c(20, 100), title='Mean Depth', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_100.dp_stats.mean), geom='line', pad=FALSE, col='grey30')
p <- p + stat_ecdf(aes(qc_50.dp_stats.mean), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.dp_stats.mean), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'dpMean_ICE_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_cumulative(df, aes(qc_150.gq_stats.mean), 'Mean Genotype Quality', 0,
    xlim=c(50, 60), title='Mean Genotype Quality', save_figure=save_figures, file="")
p <- p + stat_ecdf(aes(qc_100.gq_stats.mean), geom='line', pad=FALSE, col='grey30')
p <- p + stat_ecdf(aes(qc_50.gq_stats.mean), geom='line', pad=FALSE, col='grey50')
p <- p + stat_ecdf(aes(qc_0.gq_stats.mean), geom='line', pad=FALSE, col='grey70')
p
file=paste0(PLOTS,'gqMean_ICE_cdf')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')


df_50 <- df[,c('s', names(df)[grep("_50", names(df))])]
df_50$padding <- '50'
names(df_50) <- gsub("_50", "", names(df_50))

df_100 <- df[,c('s', names(df)[grep("_100", names(df))])]
names(df_100) <- gsub("_100", "", names(df_100))
df_100$padding <- '100'

df_0 <- df[,c('s', names(df)[grep("_0", names(df))])]
names(df_0) <- gsub("_0", "", names(df_0))
df_0$padding <- '0'

df_150 <- df[,c('s', names(df)[grep("_150", names(df))])]
names(df_150) <- gsub("_150", "", names(df_150))
df_150$padding <- '150'

df_full_ICE <- bind_rows(df_0, df_50, df_100, df_150)
df_full_ICE$padding <- ordered(df_full_ICE$padding, levels = c("0", "50", "100", "150"))

save_figures <- TRUE

create_pretty_boxplots(df_full_ICE, aes(x=padding, y=qc.n_called), aes(color=padding),
    x_label='n called', y_label='padding (bp)', title="", save_figure=save_figures,
    file=paste0(PLOTS,'n_called_ICE'), n_ticks=5)

create_pretty_boxplots(df_full_ICE, aes(x=padding, y=qc.r_ti_tv), aes(color=padding),
    x_label='Ti/Tv', y_label='padding (bp)', title="", save_figure=save_figures,
    file=paste0(PLOTS,'ti_tv_ICE'), n_ticks=5)


