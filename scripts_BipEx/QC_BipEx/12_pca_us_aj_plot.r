library(ggplot2)
library(dplyr)
library(ggsci)
library(gridExtra)

# Load plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')

# Get file locations, plotting locations, and thresholds.
source("r_functions_and_parameters/r_options_BipEx.r")
source("r_functions_and_parameters/helpful_functions.r")

save_figures <- FALSE

PCA_SCORES_USA <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_pca_scores_USA_samples.tsv"

df_PCs_USA <- fread(PCA_SCORES_USA, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df_PCs_USA <- df_PCs_USA[sample(nrow(df_PCs_USA), replace=FALSE),]

i <- 1

aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.PHENOTYPE_COARSE')
p <- create_pretty_scatter(df_PCs_USA, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1),
  n_x_ticks=5, x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.PROJECT_OR_COHORT')
create_pretty_scatter(df_PCs_USA, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1, '_batch'), n_x_ticks=5,
  x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.LOCATION')
create_pretty_scatter(df_PCs_USA, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1, '_collection'), n_x_ticks=5,
  x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))

# Now let's draw ellipses...
# Warning: Manual!
plot(df_PCs_USA$PC1, df_PCs_USA$PC2)

ellipse.1 <- ellipse(0, 0.14, -0.065, 0.02, 0.02)
ellipse.2 <- ellipse(0, 0.05, -0.04, 0.02, 0.02)

lines(ellipse.1[,1], ellipse.1[,2], lwd=2, col='indianred3')
lines(ellipse.2[,1], ellipse.2[,2], lwd=2, col='indianred3')

in_ell_1 <- in_ellipse(df_PCs_USA$PC1, df_PCs_USA$PC2, 0, 0.14, -0.065, 0.02, 0.02)
in_ell_2 <- in_ellipse(df_PCs_USA$PC1, df_PCs_USA$PC2, 0, 0.05, -0.04, 0.02, 0.02)

AJ <- 2 * in_ell_1 + in_ell_2

# Create labellings
df_PCs_USA$AJ <- factor(AJ, labels=c("not AJ", "half AJ", "AJ"))

save_figures <- TRUE
p <- create_pretty_scatter(df_PCs_USA, aes_string(x='PC1', y='PC2', color='AJ'),
	file=paste0(PLOTS, "12_AJ_in_US_PC1_PC2.jpg"),
	x_label='Principal Component 1', y_label='Principal Component 2',
	save_figure=save_figures)

m1 <- 0.55
c1 <- -0.015
m2 <- -1
c2 <- 0.185

p <- p + geom_abline(intercept=c1, slope=m1) + geom_abline(intercept=c2, slope=m2)
p

URV_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv'
PCA_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/09_pca_scores.tsv'

df_PCs <- fread(PCA_SCORES, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df_URVs <- fread(URV_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df_URVs)[names(df_URVs) == 'Samples'] <- 's'
df <- merge(df_PCs, df_URVs)
df_URVs_PCs_USA <- merge(df_PCs_USA, df, by='s', all.y=TRUE)

df_URVs_PCs_USA$AJ <- as.character(df_URVs_PCs_USA$AJ)
df_URVs_PCs_USA$AJ[is.na(df_URVs_PCs_USA$AJ)] <- 'non American'
df_URVs_PCs_USA$AJ <- factor(df_URVs_PCs_USA$AJ)

aes <- aes_string(x='PC1.y', y='n_URV_SNP', color='AJ')
p <- create_pretty_scatter(df_URVs_PCs_USA, aes, save_figure=save_figures, file='12_PC1_URVs',
  n_x_ticks=5, x_label='Principal Component 1', y_label='n URV SNPs')

df_PCs_USA$remove <- above_below_line(df_PCs_USA$PC1, df_PCs_USA$PC2, m1, c1, "below") | above_below_line(df_PCs_USA$PC1, df_PCs_USA$PC2, m2, c2, "above")

create_pretty_scatter(df_PCs_USA, aes_string(x='PC1', y='PC2', color='remove'),
	file=paste0(PLOTS, "12_AJ_remove_in_US_PC1_PC2"),
	x_label='Principal Component 1', y_label='Principal Component 2', key_label='Samples to remove',
	save_figure=save_figures)

# Add this labelling, and now include these colours on the original 1kg PC plot.
# Now, perform sanity check that these two clusters are AJs.

df_1kg <- fread(PCA_1KG_SCORES, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
    mutate(SUPER_POPULATION = factor(SUPER_POPULATION, levels=c('Control', 'Bipolar Disorder', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS')))
df_1kg <- df_1kg[sample(nrow(df_1kg), replace=FALSE),]

df_1kg <- merge(df_1kg, df_PCs_USA[c('s', 'AJ')], by='s', all=TRUE)

# Now, restrict to the strictly defined Europeans, plus these AJs - and run PCA again.
EUR_strict <- fread(EUROPEAN_SAMPLES_STRICT, sep='\t', stringsAsFactors=FALSE, header=FALSE, data.table=FALSE)$V1
EUR_strict <- setdiff(EUR_strict, (df_PCs_USA %>% filter(remove==TRUE) %>% select(s))$s)

# Now, restrict to this strict definition of the Europeans, and consider just the mainland
# Europeans and the Swedes separately.
df_out_EUR <- df_1kg %>% filter(s %in% EUR_strict) %>% select(s)
df_out_mainland_EUR <- df_1kg %>% filter(s %in% EUR_strict) %>% filter(!(LOCATION %in% c("Stockholm, SWE", "Umea, SWE"))) %>% select(s)
df_out_SWE <- df_1kg %>% filter(s %in% EUR_strict) %>% filter(LOCATION %in% c("Stockholm, SWE", "Umea, SWE")) %>% select(s)

fwrite(df_out_EUR, file=EUROPEANS, quote=FALSE, row.names=FALSE, col.names=FALSE)
fwrite(df_out_mainland_EUR, file=MAINLAND_EUROPEANS, quote=FALSE, row.names=FALSE, col.names=FALSE)
fwrite(df_out_SWE, file=SWEDES, quote=FALSE, row.names=FALSE, col.names=FALSE)

AJ <- df_PCs_USA$s[df_PCs_USA$AJ %in% c("AJ")]
df_AJ_EUR <- df_1kg[which(df_1kg$s %in% union(AJ, EUR_strict)),]
where_AJ <- which(df_AJ_EUR$s %in% AJ)
df_AJ_EUR$phenotype.PHENOTYPE_COARSE[where_AJ] <- as.character(df_AJ_EUR$AJ[where_AJ])

df_out_classify_AJ <- df_AJ_EUR %>% select(s)
fwrite(df_out_classify_AJ, file = EUROPEAN_AND_AJ_SAMPLES, quote=FALSE, row.names=FALSE, col.names=FALSE)

AJ_SAMPLES <- '../../samples_BipEx/12_aj.sample_list'
# HALF_AJ_SAMPLES <- '../../samples_BipEx/12_half_aj.sample_list'

# These can be used later - if I have time to analyse the AJs.
df_out_AJ <- data.frame(df_PCs_USA$s[df_PCs_USA$AJ == "AJ"])
fwrite(df_out_AJ, file = AJ_SAMPLES, quote=FALSE, row.names=FALSE, col.names=FALSE) 


