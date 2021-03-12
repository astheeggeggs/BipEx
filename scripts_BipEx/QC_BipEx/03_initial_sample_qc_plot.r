library(data.table)
# Plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options_BipEx.r
source("r_functions_and_parameters/r_options_BipEx.r")

QC_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv"
df <- fread(QC_FILE, stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE)

# df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_final.tsv") %>% rename(s=SAMPLE_ALIAS)

# Updated to correct psychosis and include AAO variables (April 2020)
df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv") %>% rename(s=SAMPLE_ALIAS)
df <- merge(df_pheno, df,  by='s')
# Let's rename the PI column so it's not stupid.
# Remove the 1000G samples.

# Note that 'phenotype.PHENOTYPE_COARSE etc is the old phenotype data' - or rather, we re-annotate to ensure it is completely up to date.
names(df) <- gsub("qc_padded_ice\\.", "", names(df))

# Remove the low coverage samples as defined by Excel spreadsheet e-mailed from Laura Gauthier
low_coverage <- fread("Dalio_Low_Coverage_Samples.txt")$SAMPLE

if (any(is.na(df$PI))) {
    df <- df[-which(is.na(df$PI)),]
}
df <- df[-which(df$s %in% low_coverage),]
df <- df[sample(nrow(df), replace=FALSE),]

save_figures <- TRUE
# PDFs, no splitting.
create_pretty_hist(df, aes(x=call_rate), 'Call Rate', T_sample_callRate,
    title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_hist'))
create_pretty_hist(df, aes(x=PCT_CONTAMINATION), 'Contamination', T_pct_contamination,
    binwidth=0.0001, xlim=c(0,0.1), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_hist'))
create_pretty_hist(df, aes(x=PCT_CHIMERAS), 'Chimeric Reads', T_pct_chimeras,
    binwidth=0.00002, xlim=c(0,0.02), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_hist'))
create_pretty_hist(df, aes(x=dp_stats.mean), 'Mean Depth', T_dpMean,
    binwidth=2, xlim=c(10, 150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_hist'))
create_pretty_hist(df, aes(x=gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    binwidth=0.5, xlim=c(20, 70), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_hist'))

# CDFs, no splitting.
create_pretty_cumulative(df, aes(call_rate), 'Call Rate', T_sample_callRate,
    xlim=c(0.75,1), title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_cdf'))
create_pretty_cumulative(df, aes(PCT_CONTAMINATION), 'Contamination', T_pct_contamination,
    xlim=c(0,0.1), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_cdf'))
create_pretty_cumulative(df, aes(PCT_CHIMERAS), 'Chimeric Reads', T_pct_chimeras,
    xlim=c(0,0.02), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_cdf'))
create_pretty_cumulative(df, aes(dp_stats.mean), 'Mean Depth', T_dpMean,
    xlim=c(10,150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_cdf'))
create_pretty_cumulative(df, aes(gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    xlim=c(20,70), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_cdf'))

# Split by LOCATION.
# DEV: need to be careful here - location is the location of the centre, not the location of where the sampling was done.
# For the presentation, make things larger and clearer.
legend_batch <- FALSE
legend_collection <- TRUE
legend_phenotype <- TRUE
save_figures <- TRUE
# y_label='Collection'
y_label_batch <- 'Batch'
y_label_batch <- ''
titles <- c('Call Rate',
    '% Contamination',
    '% Chimeric Reads',
    'Mean Depth (DP)',
    'Mean Genotype Quality (GQ)')
titles <- c('', '', '', '', '')
alpha <- 0.8

create_pretty_boxplots(df, aes(x=LOCATION, y=call_rate), aes(color=PROJECT_OR_COHORT),
    T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='Batch',
    xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS,'03_callRate_by_collection'), n_ticks=5, alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CONTAMINATION), aes(color=PROJECT_OR_COHORT),
    T_pct_contamination, x_label='% Contamination', y_label=y_label_batch, key_label='Batch',
    xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_batch, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS,'03_contaminiation_by_collection'), n_ticks=5, alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CHIMERAS), aes(color=PROJECT_OR_COHORT),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label=y_label_batch, key_label='Batch',
    xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_batch, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS,'03_chimeras_by_collection'), n_ticks=5, alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=dp_stats.mean), aes(color=PROJECT_OR_COHORT),
    T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='Batch',
    xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[4], save_figure=save_figures,
    file=paste0(PLOTS,'03_dpMean_by_collection'), alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=gq_stats.mean), aes(color=PROJECT_OR_COHORT),
    T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='Batch',
    xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[5], save_figure=save_figures,
    file=paste0(PLOTS,'03_gqMean_by_collection'), alpha=alpha)

create_pretty_boxplots(df, aes(x=LOCATION, y=call_rate), aes(color=PHENOTYPE_COARSE),
    T_sample_callRate, x_label='Call Rate', y_label='Collection', key_label='Phenotype',
    xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_phenotype, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS,'03_callRate_by_collection_col_CC'), alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CONTAMINATION), aes(color=PHENOTYPE_COARSE),
    T_pct_contamination, x_label='% Contamination', y_label='Collection', key_label='Phenotype',
    xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_phenotype, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS,'03_contaminiation_by_collection_col_CC'), alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CHIMERAS), aes(color=PHENOTYPE_COARSE),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label='Collection', key_label='Phenotype',
    xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_phenotype, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS,'03_chimeras_by_collection_col_CC'), alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=dp_stats.mean), aes(color=PHENOTYPE_COARSE),
    T_dpMean, x_label='Mean Depth', y_label='Collection', key_label='Phenotype',
    xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_phenotype, title=titles[4], save_figure=save_figures,
    file=paste0(PLOTS,'03_dpMean_by_collection_col_CC'), alpha=alpha)
create_pretty_boxplots(df, aes(x=LOCATION, y=gq_stats.mean), aes(color=PHENOTYPE_COARSE),
    T_gqMean, x_label='Mean Genotype Quality', y_label='Collection', key_label='Batch',
    xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_phenotype, title=titles[5], save_figure=save_figures,
    file=paste0(PLOTS,'03_gqMean_by_collection_col_CC'), alpha=alpha)

# Split by case/control
create_pretty_boxplots(df, aes(x=PHENOTYPE_COARSE, y=call_rate), aes(color=PROJECT_OR_COHORT),
    T_sample_callRate, x_label='Call Rate', y_label='Status',
    key_label='Batch', xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1],
    save_figure=save_figures, file=paste0(PLOTS,'03_callRate_by_status'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PHENOTYPE_COARSE, y=PCT_CONTAMINATION), aes(color=PROJECT_OR_COHORT),
    T_pct_contamination, x_label='% Contamination', y_label='Status',
    key_label='Batch', xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_batch, title=titles[2],
    save_figure=save_figures, file=paste0(PLOTS,'03_contamination_by_status'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PHENOTYPE_COARSE, y=PCT_CHIMERAS), aes(color=PROJECT_OR_COHORT),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label='Status',
    key_label='Batch', xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_batch, title=titles[3],
    save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_by_status'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PHENOTYPE_COARSE, y=dp_stats.mean), aes(color=PROJECT_OR_COHORT),
    T_dpMean, x_label='Mean Depth', y_label='Status',
    key_label='Batch', xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[4],
    save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_by_status'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PHENOTYPE_COARSE, y=gq_stats.mean), aes(color=PROJECT_OR_COHORT),
    T_gqMean, x_label='Mean Genotype Quality', y_label='Status',
    key_label='Batch', xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[5],
    save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_by_status'), alpha=alpha)

# Stratify by batch.
create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=call_rate), aes(color=LOCATION),
    T_sample_callRate, y_label='Batch', x_label='Call Rate', key_label='Collection',
    xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_collection, title=titles[1],
    save_figure=save_figures, file=paste0(PLOTS,'03_callRate_by_batch'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=PCT_CONTAMINATION), aes(color=LOCATION),
    T_pct_contamination, y_label='Batch', x_label='% Contamination', key_label='Collection',
    xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_collection, title=titles[2],
    save_figure=save_figures, file=paste0(PLOTS,'03_contamination_by_batch'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=PCT_CHIMERAS), aes(color=LOCATION),
    T_pct_chimeras, y_label='Batch', x_label='% Chimeric Reads', key_label='Collection',
    xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_collection, title=titles[3],
    save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_by_batch'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=dp_stats.mean), aes(color=LOCATION),
    T_dpMean, y_label='Batch', x_label='Mean Depth', key_label='Collection',
    xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_collection, title=titles[4],
    save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_by_batch'), alpha=alpha)
create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=gq_stats.mean), aes(color=LOCATION),
    T_gqMean, y_label='Batch', x_label='Mean Genotype Quality', key_label='Collection',
    xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_collection, title=titles[5],
    save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_by_batch'), alpha=alpha)

system(paste0("cp ", PLOTS, "/03_*jpg ../../site/QC_plots/sample_plots/"))
