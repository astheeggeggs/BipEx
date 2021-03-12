rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("r_functions_and_parameters/pretty_plotting.r")
source("r_functions_and_parameters/r_options_BipEx.r")

save_figure <- TRUE

SAMPLE_BEFORE_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.after.samples.tsv'

# Ensure that the most up to date phenotypes are used for the plots.

# df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_final.tsv") %>% rename(s=SAMPLE_ALIAS)

# Updated to correct psychosis and include AAO variables (April 2020)
df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv") %>% rename(s=SAMPLE_ALIAS)

if(any(duplicated(df_pheno$s))) {
    cat('remove duplicated names...\n')
    df_pheno <- df_pheno[-which(duplicated(df_pheno$s)),]
}

df_before <- fread(SAMPLE_BEFORE_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='Before Variant QC')
df_before <- merge(df_before, df_pheno, by='s')

df_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='After Variant QC')
df_after <- merge(df_after, df_pheno, by='s')

df <- bind_rows(df_before, df_after) %>%
  mutate(phase=factor(phase, levels=c('Before Variant QC', 'After Variant QC')))

## Split by Location.
# Colour by Case-status.
# Number of singletons.

y_labels <- 'Location'
y_labels <- ''
alpha <- 0.8

create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(LOCATION)),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '15_nSingletonsbyLocationColPheno'), y_label=y_labels,
    alpha=alpha)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(LOCATION)),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '15_rHetHomVarbyLocationColPheno'), y_label=y_labels,
    alpha=alpha)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(LOCATION)),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '15_rInsertionDeletionbyLocationColPheno'), y_label=y_labels,
    alpha=alpha)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(LOCATION)),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '15_rTiTvbyLocationColPheno'), n_ticks=5, y_label=y_labels,
    alpha=alpha)

## Colour by Batch.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(LOCATION)),
    aes(color=PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '15_nSingletonsbyLocationColBatch', y_label=y_labels,
        alpha=alpha))

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(LOCATION)),
    aes(color=PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '15_rHetHomVarbyLocationColBatch'), y_label=y_labels,
    alpha=alpha)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(LOCATION)),
    aes(color=PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '15_rInsertionDeletionbyLocationColBatch'), y_label=y_labels,
    alpha=alpha)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(LOCATION)),
    aes(color=PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '15_rTiTvbyLocationColBatch'), n_ticks=5, y_label=y_labels,
    alpha=alpha)

### Split by batch.
y_labels <- 'Batch'
y_labels <- ''
## Colour by Case-status.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=PROJECT_OR_COHORT),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_nSingletonsbyBatchColPheno'), y_label=y_labels,
    alpha=alpha)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=PROJECT_OR_COHORT),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rHetHomVarbyBatchColPheno'), y_label=y_labels,
    alpha=alpha)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=PROJECT_OR_COHORT),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rInsertionDeletionbyBatchColPheno'), y_label=y_labels,
    alpha=alpha)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=PROJECT_OR_COHORT),
    aes(color=PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rTiTvbyBatchColPheno'), n_ticks=5, y_label=y_labels,
    alpha=alpha)

## Colour by Location.
# Number of singletons
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_nSingletonsbyBatchColLocation'), y_label=y_labels,
    alpha=alpha)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rHetHomVarbyBatchColLocation'), y_label=y_labels,
    alpha=alpha)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rInsertionDeletionbyBatchColLocation'), y_label=y_labels,
    alpha=alpha)
# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '15_rTiTvbyBatchColLocation'), n_ticks=5, y_label=y_labels,
    alpha=alpha)

q <- function(x) {
  print(subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x))))
  y <- subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x)))
  if(length(y)==0) y <- NA
  return(y)
}

# Number of singletons
p <- create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rHetHomVar
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rInsertionDeletion
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rTiTv
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=PROJECT_OR_COHORT),
    aes(color=LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

system(paste0("cp ", PLOTS, "/15_*jpg ../../site/QC_plots/sample_plots/"))

