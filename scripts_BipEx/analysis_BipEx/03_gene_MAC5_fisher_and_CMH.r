library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggrepel)
library(rjson)

source("../QC_BipEx/r_functions_and_parameters/pretty_plotting.r")
source('r_functions/burden_tests.r')

GENE_OUT_BP_including_BPSCZ_MAC5_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool.tsv.bgz | gzcat'
GENE_OUT_SCZ_MAC5_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_SCZ_MAC5_gene_bool.tsv.bgz | gzcat'

GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool_by_location.tsv.bgz | gzcat'
GENE_OUT_SCZ_MAC5_BY_LOCATION_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_SCZ_MAC5_gene_bool_by_location.tsv.bgz | gzcat'

create_combined_qq_plot <- function(plot_name, phenotype_tests, consequences,
	dt, n_list, consequences_translation, titles, n_permutations=1000, width=4.5, height=4.5, read_dt_plot_filename=FALSE,
	Rdata_root="../../analysis_plots/gene_counts_qq/Rdata_files/", include_col=TRUE)
{
	pdf(plot_name, width=width, height=height)
	
	k <- 1
	titles <- as.vector(t(outer(titles, consequences_translation, FUN=paste, sep=":\n")))
	for (test in phenotype_tests) {
		cat(test, "...\n")
		for (consequence in consequences)
		{
			cat(consequence, "...")
			# For the Bipolar counts
			dt_plot_filename <- paste0(Rdata_root, test, "_", consequence, "_plot.Rdata")
			lookup_filename <- paste0(Rdata_root, test, "_", consequence, "_lookup.Rdata")
			
			# Obtain the number of cases and controls by summing the relevant entries of the first row.
			n_cases <- dt[1,][[paste0(test, ".case_count")]] + dt[1,][[paste0(test, ".case_no_count")]]
			n_controls <- dt[1,][[paste0(test, ".control_count")]] + dt[1,][[paste0(test, ".control_no_count")]]
			
			if (!(read_dt_plot_filename & file.exists(dt_plot_filename))) {
				dt_list <- run_qq_create_fisher_and_perm(
					dt, n_cases, n_controls,
					paste0(test, ".case_count"), paste0(test, ".control_count"),
					n_permutations=n_permutations, consequence=consequence,
					save_lookup=TRUE, lookup_filename=lookup_filename,
					read_lookup=FALSE, include_col=include_col
				)
				dt_plot <- dt_list$dt_plot
				results_dt <- dt_list$results_dt
				save(dt_plot, results_dt, file=dt_plot_filename)
			}
			create_qq_plot_and_table(dt_plot_filename, qq_labels=20, plot_title=titles[k],
				table_out=paste0("../../analysis_plots/gene_counts_qq/plots/", test))
			k <- k+1
		}
		cat("\n")
	}

	dev.off()
}

create_fisher_supplementary_table <- function(dt, phenotype_tests, log=FALSE)
{
	init <- TRUE
	for (phenotype in phenotype_tests)
	{	
		cat(paste0(phenotype, "..."))
		cat("obtaining observed p-values...\n")
		dt_obs <- get_fisher_genes(phenotype, dt) %>% 
		select(gene_symbol, consequence_category, pval, OR)

		if (log == TRUE) {
			dt_obs <- dt_obs %>% mutate(pval = log10(pval), log_OR = log10(OR)) %>% select(-OR)
			names_to_append <- which(names(dt_obs) %in% c("pval", "log_OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("log_pval", "log_OR"))
		} else {
			names_to_append <- which(names(dt_obs) %in% c("pval", "OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("pval", "OR"))
		}
		cat("obtained observed p-values...\n")
		
		if (init) {
			dt_full <- merge(dt, dt_obs, by=c("gene_symbol", "consequence_category"))
			init <- FALSE
		} else {
			dt_full <- merge(dt_full, dt_obs, by=c("gene_symbol", "consequence_category"))
		}
	}
	return(data.table(dt_full))
}

create_CMH_supplementary_table <- function(dt, phenotype_tests, log=FALSE)
{
	init <- TRUE
	for (phenotype in phenotype_tests)
	{	
		cat(paste0(phenotype, "..."))
		cat("obtaining observed p-values...\n")
		dt_obs <- get_CMH_genes(phenotype, dt) %>% 
		select(gene_symbol, consequence_category, pval_cc, R)

		if (log == TRUE) {
			dt_obs <- dt_obs %>% 
				mutate(pval = -log10(pval_cc), log_OR = log10(R)) %>% 
				select(-R, -pval_cc)
			names_to_append <- which(names(dt_obs) %in% c("pval", "log_OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("log_pval", "log_OR"))
		} else {
			dt_obs <- dt_obs %>% rename(pval = pval_cc, OR = R)
			names_to_append <- which(names(dt_obs) %in% c("pval", "OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("pval", "OR"))
		}
		cat("obtained observed p-values...\n")
		
		if (init) {
			dt_full <- dt_obs
			init <- FALSE
		} else {
			dt_full <- merge(dt_full, dt_obs, by=c("gene_symbol", "consequence_category"))
		}
	}
	return(data.table(dt_full))
}

create_combined_qq_plots <- function(dt, consequences, titles, consequences_translation,
	phenotype_tests, n_perms=10, qq_labels=20, plot_name='plot.pdf', width=4.5, height=4,
	include_col=TRUE, separate_plots=FALSE, separate_plot_name=NULL, title.hjust=0)
{
	k <- 1
	titles <- as.vector(t(outer(titles, consequences_translation, FUN=paste, sep=":\n")))
	pdf(plot_name, width=width, height=height)

	for (phenotype in phenotype_tests)
	{	
		cat(paste0(phenotype, "..."))
		cat("obtaining observed p-values...\n")
		dt_obs <- get_fisher_genes(phenotype, dt) %>% rename(labels=gene_symbol) %>% 
		select(labels, consequence_category, pval, OR) %>% 
		mutate(
			pval = -log10(pval),
			log_OR = log10(OR)
			)
		cat("obtained observed p-values...\n")

		dt_exp <- create_permutation_genes(n_perms, phenotype, dt) %>% 
		mutate(p_perm = -log10(pval))

		dt_plot <- data.table(dt_obs)[!is.na(pval),]
		dt_exp <- data.table(dt_exp)[!is.na(p_perm)]

		setkeyv(dt_plot, c("consequence_category", "pval"))
		setkeyv(dt_exp, c("consequence_category", "p_perm"))
		dt_plot[, p_perm := dt_exp$p_perm]

		cat("Creating QQ plot...\n")

		save(dt_plot, file=paste0(gsub(".pdf", "", plot_name), '_', phenotype, '.Rdata'))

		for (consequence in consequences)
		{	
			dt_plot_current <- dt_plot[consequence_category == consequence, ]

			if (include_col)
			{
				max_replace <- max(abs(dt_plot_current$log_OR[is.finite(dt_plot_current$log_OR)]))
				dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)] <- sign(dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)]) * max_replace

				p <- create_pretty_qq_plot(
					dt_plot_current, aes(x=p_perm, y=pval, color=log_OR),
					n_to_include=qq_labels, cex_label=2, plot_title=titles[k], gradient=include_col, 
					gradient_title=TeX("$\\log(OR)$"), save_figure=separate_plots,
					file=paste0(separate_plot_name, '_', phenotype, '_', consequence),
					width=width*2.54*10, height=height*2.54*10, title.hjust=title.hjust)
			} else {
				p <- create_pretty_qq_plot(
					dt_plot_current, aes(x=p_perm, y=pval),
					n_to_include=qq_labels, cex_label=2, plot_title=titles[k], gradient=include_col,
					save_figure=separate_plots,
					file=paste0(separate_plot_name, '_', phenotype, '_', consequence),
					width=width*2.54*10, height=height*2.54*10, title.hjust=title.hjust)
			}
			k <- k+1
		}
		cat("Created QQ plot...\n")
	}
	dev.off()
	return(dt_plot)
}

create_combined_CMH_qq_plots <- function(dt, consequences, titles, consequences_translation,
	phenotype_tests, n_perms=10, qq_labels=20, plot_name='plot.pdf', width=4.5, height=4,
	include_col=TRUE, separate_plots=FALSE, separate_plot_name=NULL, title.hjust=0, force=TRUE)
{
	k <- 1
	titles <- as.vector(t(outer(titles, consequences_translation, FUN=paste, sep=":\n")))
	pdf(plot_name, width=width, height=height)
	
	for (phenotype in phenotype_tests)
	{
		cat(paste0(phenotype, "..."))
		
		Rdata_file <- paste0(gsub(".pdf", "", plot_name), '_', phenotype, '.Rdata')
		if (file.exists(Rdata_file) & !force)
		{
			cat(paste0("Warning: using existing Rdata file to create plots. If you do not want to use this existing",
				" data, remove the file located at ", Rdata_file, "\n"))
			load(Rdata_file)

		} else {
			cat("obtaining observed p-values...\n")
			dt_obs <- get_CMH_genes(phenotype, dt) %>% rename(labels=gene_symbol) %>% 
			select(labels, consequence_category, pval_cc, R) %>% 
			mutate(
				pval = -log10(pval_cc),
				log_OR = log10(R)
				)
			cat("obtained observed p-values...\n")
			dt_exp <- create_permutation_CMH_genes(n_perms, phenotype, dt) %>% 
			mutate(p_perm = -log10(pval_cc))

			# Number that are undefined should be the same in the permuations and real data
			# these are the genes for which there are no case or control counts.
			dt_plot <- data.table(dt_obs)[!is.na(pval_cc),]
			dt_exp <- data.table(dt_exp)[!is.na(p_perm)]
			setkeyv(dt_plot, c("consequence_category", "pval"))
			setkeyv(dt_exp, c("consequence_category", "p_perm"))
			dt_plot[, p_perm := dt_exp$p_perm]

			cat("Creating QQ plot...\n")
			save(dt_plot, file=paste0(gsub(".pdf", "", plot_name), '_', phenotype, '.Rdata'))
			cat("Saved .Rdata file\n")
		}

		for (consequence in consequences)
		{	
			dt_plot_current <- dt_plot[consequence_category == consequence, ]

			if (include_col) {

				max_replace <- max(abs(dt_plot_current$log_OR[is.finite(dt_plot_current$log_OR)]))
				dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)] <- sign(dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)]) * max_replace

				p <- create_pretty_qq_plot(
					dt_plot_current, aes(x=p_perm, y=pval, color=log_OR),
					n_to_include=qq_labels, cex_label=2, plot_title=titles[k], gradient=include_col, 
					gradient_title=TeX("$\\log(OR)$"), save_figure=separate_plots,
					file=paste0(separate_plot_name, '_', phenotype, '_', consequence),
					width=width*2.54*10, height=height*2.54*10, title.hjust=title.hjust)
			} else {
				p <- create_pretty_qq_plot(
					dt_plot_current, aes(x=p_perm, y=pval),
					n_to_include=qq_labels, cex_label=2, plot_title=titles[k], gradient=include_col,
					save_figure=separate_plots,
					file=paste0(separate_plot_name, '_', phenotype, '_', consequence),
					width=width*2.54*10, height=height*2.54*10, title.hjust=title.hjust)
			}
			k <- k+1
		}
		cat("Created QQ plot...\n")
	}
	dev.off()
}

website_qq_plots <- function(Rdata, consequences, consequences_translation,
	plot_name, qq_labels=20, width=4.5, height=4, include_col=TRUE, 
	title.hjust=0, include_title=TRUE, include_legend=TRUE, legend_max_min=NULL)
{
	k <- 1
	if (include_title) {
		titles <- consequences_translation
	} else {
		titles <- rep("", length(consequences_translation))
	}

	load(Rdata)

	cat("Creating QQ plot...\n")
	for (consequence in consequences)
	{	
		cat("consequence:", consequence, "\n")
		dt_plot_current <- dt_plot[consequence_category == consequence, ]

		if (include_col) {

			max_replace <- max(abs(dt_plot_current$log_OR[is.finite(dt_plot_current$log_OR)]))
			dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)] <- sign(dt_plot_current$log_OR[!is.finite(dt_plot_current$log_OR)]) * max_replace

			p <- create_pretty_qq_plot(
				dt_plot_current, aes(x=p_perm, y=pval, color=log_OR),
				n_to_include=qq_labels, plot_title=titles[k], gradient=include_col, 
				gradient_title=TeX("$\\log(OR)$"), save_figure=FALSE, print_p=FALSE, title.hjust=title.hjust,
				cex_labels=2.5)
			if (include_legend == FALSE) {
				p <- p + theme(legend.position="none")
			}
			ggsave(paste0(plot_name, '_', consequence, '.jpg'), p, width=width, height=height)

		} else {
			p <- create_pretty_qq_plot(
				dt_plot_current, aes(x=p_perm, y=pval),
				n_to_include=qq_labels, plot_title=titles[k], gradient=include_col,
				save_figure=FALSE, print_p=FALSE, file=paste0(plot_name, '_', consequence), title.hjust=title.hjust,
				cex_labels=2.5)
			ggsave(paste0(plot_name, '_', consequence, '.jpg'), p, width=width, height=height)
		}
		k <- k+1
	}
	cat("Created QQ plot...\n")

}

wrapper_for_figures <- function(
	qq_function,
	dt_location,
	dir_figures_output, figures_output,
	dir_separate_figures_output, separate_figures_output,
	phenotype_tests, titles, consequences, consequences_translation,
	n_perms=10, qq_labels=20, width=4.5, height=4, include_col=TRUE,
	separate_plots=TRUE) 
{
	if ((!file.exists(figures_output)) | (!file.exists(dir_separate_figures_output)))
	{
		dt <- fread(dt_location)
		names(dt) <- gsub("burden.", "", names(dt)) 
		# The last row contains only NAs
		dt <- dt[-nrow(dt),]
		dir.create(dir_figures_output, showWarnings = FALSE)
		dir.create(dir_separate_figures_output, showWarnings=FALSE)
		qq_function(dt, consequences, titles, consequences_translation,
			phenotype_tests, n_perms=n_perms, qq_labels=qq_labels, plot_name=figures_output,
			width=width, height=height, include_col=TRUE, separate_plots=separate_plots,
			separate_plot_name=separate_figures_output)
	}

	# Write the table
	# tables_output <- gsub(".pdf", ".tsv", figures_output)
	# if (!file.exists(tables_output)) {
	# 	dt <- fread(dt_location)
	# 	names(dt) <- gsub("burden.", "", names(dt)) 
	# 	# The last row contains only NAs
	# 	dt <- dt[-nrow(dt),]
	# 	dt_full <- create_fisher_supplementary_table(dt, phenotype_tests)
	# 	fwrite(dt_full, file=tables_output, sep="\t")
	# }
}

consequences <- c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv")
consequences_translation <- c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV")

## Bipolar Disorder

BP_phenotype_tests <- c("is_BP1", "is_BP2", "is_BP", "is_BPPSY", "is_BP_no_PSY")
BP_gnom_non_psych_phenotype_tests <- paste0("gnom_non_psych.", BP_phenotype_tests)
BP_phenotype_tests_list <- list(BP_phenotype_tests, BP_gnom_non_psych_phenotype_tests)

BP_titles <- c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder", "Bipolar Disorder including Schizoaffective",
	"Bipolar Disorder NOS", "Schizoaffective", "Bipolar Disorder with psychosis", "Bipolar Disorder without psychosis")

BP_dt <- c(GENE_OUT_BP_including_BPSCZ_TSV, GENE_OUT_BP_including_BPSCZ_TSV,
	GENE_OUT_BP_including_BPSCZ_MAC5_TSV, GENE_OUT_BP_including_BPSCZ_MAC5_TSV)

BP_dir_figures <- c("../../analysis_plots/gene_counts_qq/plots", "../../analysis_plots/gene_counts_qq/plots")
BP_plot_names <- paste0(
	BP_dir_figures,
	c("/BP_gene_fisher_MAC5_qq.pdf", "/BP_gene_fisher_MAC5_gnom_non_psych_qq.pdf")
	)

BP_dir_separate_figures_output <- c(
	"../../analysis_plots/gene_counts_qq/plots/fisher_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/fisher_gnom_non_psych_MAC5"
	)

BP_separate_figures_output <- paste0(
	BP_dir_separate_figures_output,
	c("/BP_gene_fisher_MAC5_qq", "/BP_gene_fisher_gnom_non_psych_MAC5_qq")
	)

## Schizophrenia

SCZ_phenotype_tests <- c("is_SCZ")
SCZ_gnom_non_psych_phenotype_tests <- paste0("gnom_non_psych.", SCZ_phenotype_tests)
SCZ_phenotype_tests_list <- list(SCZ_phenotype_tests, SCZ_gnom_non_psych_phenotype_tests)

SCZ_titles <- c("Schizophrenia")
SCZ_dt <- c(GENE_OUT_SCZ_MAC5_TSV, GENE_OUT_SCZ_MAC5_TSV)
SCZ_dir_figures <- BP_dir_figures

SCZ_plot_names <- paste0(
	SCZ_dir_figures,
	c("/SCZ_gene_fisher_MAC5_qq.pdf", "/SCZ_gene_fisher_MAC5_gnom_non_psych_qq.pdf")
	)

SCZ_dir_separate_figures_output <- BP_dir_separate_figures_output

SCZ_separate_figures_output <- paste0(
	SCZ_dir_separate_figures_output,
	c("/SCZ_gene_fisher_MAC5_qq", "/SCZ_gene_fisher_gnom_non_psych_MAC5_qq")
	)

## Bipolar Disorder
for (i in 1:length(BP_plot_names)) {
	wrapper_for_figures(create_combined_qq_plots, BP_dt[i], BP_dir_figures[i], BP_plot_names[i],
		BP_dir_separate_figures_output[i], BP_separate_figures_output[i], BP_phenotype_tests_list[[i]],
		BP_titles, consequences, consequences_translation, n_perms=20)
}

## Schizophrenia
for (i in 1:length(SCZ_plot_names)) {
	wrapper_for_figures(create_combined_qq_plots, SCZ_dt[i], SCZ_dir_figures[i], SCZ_plot_names[i],
		SCZ_dir_separate_figures_output[i], SCZ_separate_figures_output[i], SCZ_phenotype_tests_list[[i]],
		SCZ_titles, consequences, consequences_translation, n_perms=20)
}

# For the website (MAC 5)
Rdatas <- paste0(
	BP_dir_figures[1],
	c(
		"/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP.Rdata",
		"/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP1.Rdata",
		"/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP2.Rdata",
		"/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BPPSY.Rdata",
		"/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP_no_PSY.Rdata"
	)
)

website_plots <- "../../analysis_plots/gene_counts_qq/website_plots/"
dir.create(website_plots, showWarnings = FALSE)
plot_names <- paste0(
	website_plots,
	c(
	  "fisher_BP_MAC5_gnom_non_psych",
	  "fisher_BP1_MAC5_gnom_non_psych",
	  "fisher_BP2_MAC5_gnom_non_psych",
	  "fisher_BPPSY_MAC5_gnom_non_psych",
	  "fisher_BP_no_PSY_MAC5_gnom_non_psych"
	)
)

include_title <- rep(c(TRUE, FALSE), 8)
include_legend <- rep(c(FALSE, TRUE), 8)
width <- rep(c(3.7, 4.5), 8)

for (i in 1:length(Rdatas)) {
	website_qq_plots(Rdatas[i], consequences, consequences_translation,
		plot_names[i], qq_labels=20, width=width[i], height=4, include_col=TRUE, 
		title.hjust=0, include_title=include_title[i], include_legend=include_legend[i],
		legend_max_min=NULL)
}

# CMH
## Bipolar Disorder

BP_dt <- c(GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV, GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV)

BP_phenotype_tests <- c("is_BP1", "is_BP2", "is_BP", "is_BPPSY", "is_BP_no_PSY")
BP_gnom_non_psych_phenotype_tests <- paste0("gnom_non_psych.", BP_phenotype_tests)
BP_phenotype_tests_list <- list(
	BP_phenotype_tests,
	BP_gnom_non_psych_phenotype_tests
	)

BP_titles <- c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder",
	"Bipolar Disorder with psychosis", "Bipolar Disorder without psychosis")
BP_plot_names <- paste0(
	BP_dir_figures,
	c("/BP_gene_CMH_MAC5_qq.pdf", "/BP_gene_CMH_MAC5_gnom_non_psych_qq.pdf")
	)
BP_dir_separate_figures_output <- c(
	"../../analysis_plots/gene_counts_qq/plots/CMH_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/CMH_gnom_non_psych_MAC5"
	)
BP_separate_figures_output <- paste0(
	BP_dir_separate_figures_output,
	c("/BP_gene_CMH_MAC5_qq", "/BP_gene_CMH_gnom_non_psych_MAC5_qq")
	)

## Schizophrenia

SCZ_dt <- c(GENE_OUT_SCZ_MAC5_BY_LOCATION_TSV, GENE_OUT_SCZ_MAC5_BY_LOCATION_TSV)
SCZ_dir_figures <- BP_dir_figures
SCZ_plot_names <- paste0(SCZ_dir_figures,
	c("/SCZ_gene_CMH_MAC5_qq.pdf", "/SCZ_gene_CMH_MAC5_gnom_non_psych_qq.pdf")
	)
SCZ_dir_separate_figures_output <- c(
	"../../analysis_plots/gene_counts_qq/plots/CMH_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/CMH_gnom_non_psych_MAC5"
)
SCZ_separate_figures_output <- paste0(
	SCZ_dir_separate_figures_output,
	c("/SCZ_gene_CMH_MAC5_qq", "/SCZ_gene_CMH_gnom_non_psych_MAC5_qq")
	)

## Bipolar Disorder
for (i in 1:length(BP_plot_names)) {
	wrapper_for_figures(create_combined_CMH_qq_plots, BP_dt[i], BP_dir_figures[i], BP_plot_names[i],
		BP_dir_separate_figures_output[i], BP_separate_figures_output[i], BP_phenotype_tests_list[[i]],
		BP_titles, consequences, consequences_translation, n_perms=20)
}

## Schizophrenia
for (i in 1:length(SCZ_plot_names)) {
	wrapper_for_figures(create_combined_CMH_qq_plots, SCZ_dt[i], SCZ_dir_figures[i], SCZ_plot_names[i],
		SCZ_dir_separate_figures_output[i], SCZ_separate_figures_output[i], SCZ_phenotype_tests_list[[i]],
		SCZ_titles, consequences, consequences_translation, n_perms=20)
}

# For the website (MAC <= 5)
Rdatas <- paste0(
	BP_dir_figures[1],
	c(
		"/BP_gene_CMH_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP.Rdata",
		"/BP_gene_CMH_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP1.Rdata",
		"/BP_gene_CMH_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP2.Rdata",
		"/BP_gene_CMH_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BPPSY.Rdata",
		"/BP_gene_CMH_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP_no_PSY.Rdata"
	)
)

website_plots <- "../../analysis_plots/gene_counts_qq/website_plots/"
dir.create(website_plots, showWarnings = FALSE)
plot_names <- paste0(
	website_plots,
	c(
	  "CMH_BP_MAC5_gnom_non_psych",
	  "CMH_BP1_MAC5_gnom_non_psych",
	  "CMH_BP2_MAC5_gnom_non_psych",
	  "CMH_BPPSY_MAC5_gnom_non_psych",
	  "CMH_BP_no_PSY_MAC5_gnom_non_psych"
	)
)

include_title <- rep(c(TRUE, FALSE), 6)
include_legend <- rep(c(FALSE, TRUE), 6)
width <- rep(c(3.7, 4.5), 6)

for (i in 1:length(Rdatas)) {
	website_qq_plots(Rdatas[i], consequences, consequences_translation,
		plot_names[i], qq_labels=20, width=width[i], height=4, include_col=TRUE, 
		title.hjust=0, include_title=include_title[i], include_legend=include_legend[i],
		legend_max_min=NULL)
}

# For the website (singletons)
Rdatas <- paste0(
	BP_dir_figures[1],
	c(
		"/BP_gene_CMH_gnom_non_psych_qq_gnom_non_psych.is_BP.Rdata",
		"/BP_gene_CMH_gnom_non_psych_qq_gnom_non_psych.is_BP1.Rdata",
		"/BP_gene_CMH_gnom_non_psych_qq_gnom_non_psych.is_BP2.Rdata",
		"/BP_gene_CMH_gnom_non_psych_qq_gnom_non_psych.is_BPPSY.Rdata",
		"/BP_gene_CMH_gnom_non_psych_qq_gnom_non_psych.is_BP_no_PSY.Rdata"
	)
)

website_plots <- "../../analysis_plots/gene_counts_qq/website_plots/"
dir.create(website_plots, showWarnings = FALSE)
plot_names <- paste0(
	website_plots,
	c(
	  "CMH_BP_gnom_non_psych",
	  "CMH_BP1_gnom_non_psych",
	  "CMH_BP2_gnom_non_psych",
	  "CMH_BPPSY_gnom_non_psych",
	  "CMH_BP_no_PSY_gnom_non_psych"
	)
)

include_title <- rep(c(TRUE, FALSE), 6)
include_legend <- rep(c(FALSE, TRUE), 6)
width <- rep(c(3.7, 4.5), 6)

for (i in 1:length(Rdatas)) {
	website_qq_plots(Rdatas[i], consequences, consequences_translation,
		plot_names[i], qq_labels=20, width=width[i], height=4, include_col=TRUE, 
		title.hjust=0, include_title=include_title[i], include_legend=include_legend[i],
		legend_max_min=NULL)
}

# Large genes to test the null

consequences <- "synonymous"
consequences_translation <- "Synonymous"

# Fisher
## Bipolar Disorder
BP_phenotype_tests <- c("is_BP")
BP_gnom_non_psych_phenotype_tests <- paste0("gnom_non_psych.", BP_phenotype_tests)
BP_phenotype_tests_list <- list(
	BP_phenotype_tests, BP_gnom_non_psych_phenotype_tests,
	BP_phenotype_tests, BP_gnom_non_psych_phenotype_tests)
BP_titles <- c("Bipolar Disorder")
BP_dir_figures <- "../../analysis_plots/gene_counts_qq/plots/null_test"
BP_plot_names <- paste0(BP_dir_figures,
	c("/BP_gene_fisher_MAC5_qq_50.pdf", "/BP_gene_fisher_MAC5_gnom_non_psych_qq_50.pdf",
	  "/BP_gene_fisher_MAC5_qq_20.pdf", "/BP_gene_fisher_MAC5_gnom_non_psych_qq_20.pdf"))
BP_dir_separate_figures_output <- c(
	"../../analysis_plots/gene_counts_qq/plots/null_test/fisher_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/fisher_gnom_non_psych_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/fisher_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/fisher_gnom_non_psych_MAC5")
BP_separate_figures_output <- paste0(BP_dir_separate_figures_output,
	c("/BP_gene_fisher_MAC5_qq_50","/BP_gene_fisher_gnom_non_psych_MAC5_qq_50",
	  "/BP_gene_fisher_MAC5_qq_20","/BP_gene_fisher_gnom_non_psych_MAC5_qq_20"))

# First, read in and restrict to the collection of genes that have at least 50/20 individuals with a MAC 5 variant in that gene.
# Read in, filter and write.
dt <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_TSV)
dt <- dt %>% filter(consequence_category == "synonymous") %>% filter((burden_gnom_non_psych.is_BP.case_count + burden_gnom_non_psych.is_BP.control_count) > 20)
fwrite(dt, file="null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous.tsv" , sep='\t')
dt <- dt %>% filter(consequence_category == "synonymous") %>% filter((burden_gnom_non_psych.is_BP.case_count + burden_gnom_non_psych.is_BP.control_count) > 50)
fwrite(dt, file="null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous.tsv" , sep='\t')
BP_dt <- c("null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous.tsv", "null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous.tsv",
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous.tsv", "null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous.tsv")

for (i in 1:length(BP_plot_names)) {
	wrapper_for_figures(create_combined_qq_plots, BP_dt[i], BP_dir_figures[i], BP_plot_names[i],
		BP_dir_separate_figures_output[i], BP_separate_figures_output[i], BP_phenotype_tests_list[[i]],
		BP_titles, consequences, consequences_translation, n_perms=20)
}

# And do the same thing for CMH.
# CMH
## Bipolar Disorder

BP_phenotype_tests <- c("is_BP")
BP_gnom_non_psych_phenotype_tests <- paste0("gnom_non_psych.", BP_phenotype_tests)
BP_phenotype_tests_list <- list(
	BP_phenotype_tests, BP_gnom_non_psych_phenotype_tests,
	BP_phenotype_tests, BP_gnom_non_psych_phenotype_tests)
BP_titles <- c("Bipolar Disorder")
BP_dir_figures <- "../../analysis_plots/gene_counts_qq/plots/null_test"
BP_plot_names <- paste0(BP_dir_figures,
	c("/BP_gene_CMH_MAC5_qq_50.pdf", "/BP_gene_CMH_MAC5_gnom_non_psych_qq_50.pdf",
	  "/BP_gene_CMH_MAC5_qq_20.pdf", "/BP_gene_CMH_MAC5_gnom_non_psych_qq_20.pdf"))
BP_dir_separate_figures_output <- c(
	"../../analysis_plots/gene_counts_qq/plots/null_test/CMH_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/CMH_gnom_non_psych_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/CMH_MAC5",
	"../../analysis_plots/gene_counts_qq/plots/null_test/CMH_gnom_non_psych_MAC5")
BP_separate_figures_output <- paste0(BP_dir_separate_figures_output,
	c("/BP_gene_CMH_MAC5_qq_50", "/BP_gene_CMH_gnom_non_psych_MAC5_qq_50",
	  "/BP_gene_CMH_MAC5_qq_20", "/BP_gene_CMH_gnom_non_psych_MAC5_qq_20"))

dt <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV) %>% filter(consequence_category == "synonymous")
genes_to_include <- dt %>% group_by(gene_symbol) %>% summarise(n_AC=sum(burden_gnom_non_psych.is_BP.case_count + burden_gnom_non_psych.is_BP.control_count))
genes_to_include_50 <- (genes_to_include %>% filter(n_AC > 50))$gene_symbol
genes_to_include_20 <- (genes_to_include %>% filter(n_AC > 20))$gene_symbol

dt <- dt %>% filter(consequence_category == "synonymous") %>% filter(gene_symbol %in% genes_to_include_20)
fwrite(dt, file="null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous_by_location.tsv" , sep='\t')
dt <- dt %>% filter(consequence_category == "synonymous") %>% filter(gene_symbol %in% genes_to_include_50)
fwrite(dt, file="null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous_by_location.tsv" , sep='\t')
BP_dt <- c("null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous_by_location.tsv", "null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous_by_location.tsv",
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous_by_location.tsv", "null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous_by_location.tsv")

## Bipolar Disorder
for (i in 1:length(BP_plot_names)) {
	wrapper_for_figures(create_combined_CMH_qq_plots, BP_dt[i], BP_dir_figures[i], BP_plot_names[i],
		BP_dir_separate_figures_output[i], BP_separate_figures_output[i], BP_phenotype_tests_list[[i]],
		BP_titles, consequences, consequences_translation, n_perms=20)
}


