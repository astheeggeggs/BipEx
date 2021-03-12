rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("r_functions_and_parameters/pretty_plotting.r")
source("r_functions_and_parameters/r_options_BipEx.r")

save_figure <- TRUE

SAMPLE_QC_IN_TARGET = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/18_final_qc.after_in_target.samples.tsv'

df <- fread(SAMPLE_QC_IN_TARGET, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df) <- gsub("imputesex.", "", names(df))
df_before <- df[,-c(grep("sample_qc_in_target", names(df)))] %>% mutate(phase="Padded target intervals")
df_after <-  df[,-c(grep("sample_qc\\.", names(df)))] %>% mutate(phase='Target intervals')
names(df_after) <- gsub("_in_target", "", names(df_after))
df <- bind_rows(df_before, df_after) %>%
  mutate(phase=factor(phase, levels=c('Padded target intervals', 'Target intervals')))

## Split by Location.
# Colour by Case-status.
# Number of singletons.

y_labels <- 'Location'
y_labels <- ''

create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyLocationColPheno'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyLocationColPheno'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyLocationColPheno'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyLocationColPheno'), n_ticks=5, y_label=y_labels)

## Colour by Batch.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyLocationColBatch', y_label=y_labels))

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyLocationColBatch'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyLocationColBatch'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(phenotype.LOCATION)),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyLocationColBatch'), n_ticks=5, y_label=y_labels)

### Split by batch.
y_labels <- 'Batch'
y_labels <- ''
## Colour by Case-status.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyBatchColPheno'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyBatchColPheno'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyBatchColPheno'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.PHENOTYPE_COARSE), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyBatchColPheno'), n_ticks=5, y_label=y_labels)

## Colour by Location.
# Number of singletons
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyBatchColLocation'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyBatchColLocation'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyBatchColLocation'), y_label=y_labels)
# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.PROJECT_OR_COHORT),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyBatchColLocation'), n_ticks=5, y_label=y_labels)

### Split by case-status.
y_labels <- 'Case-status'
y_labels <- ''
## Colour by Case-status.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyStatusColPheno'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyStatusColPheno'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyStatusColPheno'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.PROJECT_OR_COHORT), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyStatusColPheno'), n_ticks=5, y_label=y_labels)

## Colour by Location.
# Number of singletons
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_nSingletonsbyStatusColLocation'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rHetHomVarbyStatusColLocation'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rInsertionDeletionbyStatusColLocation'), y_label=y_labels)
# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.PHENOTYPE_COARSE),
    aes(color=phenotype.LOCATION), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '18_in_target_rTiTvbyStatusColLocation'), n_ticks=5, y_label=y_labels)




