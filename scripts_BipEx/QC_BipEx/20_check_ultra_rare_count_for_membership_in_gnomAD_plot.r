library(hexbin)
library(gridExtra)
library(ggplot2)
library(ggExtra)
library(data.table)

# Load plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')

# Get file locations and plotting locations.
source("r_functions_and_parameters/r_options_BipEx.r")

NOT_IN_GNOMAD_URV_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/20_URVs_not_in_gnomAD.tsv'
NOT_IN_GNOMAD_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/20_not_in_gnomAD.tsv'

df <- fread(NOT_IN_GNOMAD_URV_FILE, sep='\t', header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df), replace=FALSE),]

# Scatters of URV-SNPs against URV-Indels.
d <- ggplot(df, aes(x=n_URV_SNP, y=n_URV_indel, colour=phenotype.PROJECT_OR_COHORT)) + 
geom_point(size=0.5, alpha=0.5) + 
scale_color_d3('category10') + theme_minimal() + labs(x='n URV SNPs', y='n URV indels', color='Batch')
d <- ggExtra::ggMarginal(d, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))
print(d)

d <- ggplot(df, aes(x=n_URV_SNP, y=n_URV_indel, colour=phenotype.PROJECT_OR_COHORT)) + 
geom_point(size=0.5, alpha=0.5) + scale_color_d3('category10') + theme_minimal() + xlim(c(0,700)) + ylim(c(0,100)) +
labs(x='n URV SNPs', y='n URV indels', color='Batch')
d <- ggExtra::ggMarginal(d, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))

save_figures <- TRUE

y_labels <- c('Batch', 'Location')
y_label_batch <- c('', '')
titles <- c('Number of Singletons split by Batch and coloured by Location',
	'Number of Singletons split by Location and coloured by Phenotype')
titles <- c('', '')

create_pretty_boxplots(df, aes(x=phenotype.PROJECT_OR_COHORT, y=n_URV_SNP+n_URV_indel), aes(colour=factor(phenotype.LOCATION)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_URVs_by_batch'))
create_pretty_boxplots(df, aes(x=phenotype.LOCATION, y=n_URV_SNP+n_URV_indel), aes(colour=factor(phenotype.PROJECT_OR_COHORT)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_URVs_by_location'))

# Actually want to plot the PTV stuff - these should certainly be in the coding regions.
create_pretty_boxplots(df, aes(x=phenotype.PROJECT_OR_COHORT, y=n_URV_PTV), aes(colour=factor(phenotype.LOCATION)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_URV_PTVs_by_batch'))
create_pretty_boxplots(df, aes(x=phenotype.LOCATION, y=n_URV_PTV), aes(colour=factor(phenotype.PROJECT_OR_COHORT)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_URV_PTVs_by_location'))

df <- fread(NOT_IN_GNOMAD_FILE, sep='\t', header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df), replace=FALSE),]

# Scatters of URV-SNPs against URV-Indels.
d <- ggplot(df, aes(x=n_SNP, y=n_indel, colour=phenotype.PROJECT_OR_COHORT)) + 
geom_point(size=0.5, alpha=0.5) + 
scale_color_d3('category10') + theme_minimal() + labs(x='n URV SNPs', y='n URV indels', color='Batch')
d <- ggExtra::ggMarginal(d, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))
print(d)

d <- ggplot(df, aes(x=n_SNP, y=n_indel, colour=phenotype.PROJECT_OR_COHORT)) + 
geom_point(size=0.5, alpha=0.5) + scale_color_d3('category10') + theme_minimal() + xlim(c(0,700)) + ylim(c(0,100)) +
labs(x='n SNPs', y='n indels', color='Batch')
d <- ggExtra::ggMarginal(d, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))

save_figures <- TRUE

y_labels <- c('Batch', 'Location')
y_label_batch <- c('', '')
titles <- c('Number of Singletons split by Batch and coloured by Location',
	'Number of Singletons split by Location and coloured by Phenotype')
titles <- c('', '')

create_pretty_boxplots(df, aes(x=phenotype.PROJECT_OR_COHORT, y=n_SNP+n_indel), aes(colour=factor(phenotype.LOCATION)), 
	y_label=y_label_batch[1], x_label='Number of variants not in gnomAD', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_by_batch'))
create_pretty_boxplots(df, aes(x=phenotype.LOCATION, y=n_SNP+n_indel), aes(colour=factor(phenotype.PROJECT_OR_COHORT)),
	y_label=y_label_batch[2], x_label='Number of variants not in gnomAD', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_by_location'))

# Actually want to plot the PTV stuff - these should certainly be in the coding regions.
create_pretty_boxplots(df, aes(x=phenotype.PROJECT_OR_COHORT, y=n_PTV), aes(colour=factor(phenotype.LOCATION)), 
	y_label=y_label_batch[1], x_label='Number of PTVs not in gnomAD', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_PTVs_by_batch'))
create_pretty_boxplots(df, aes(x=phenotype.LOCATION, y=n_PTV), aes(colour=factor(phenotype.PROJECT_OR_COHORT)),
	y_label=y_label_batch[2], x_label='Number of PTVs not in gnomAD', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'20_PTVs_by_location'))

# Finally, should further split by case/control status.

