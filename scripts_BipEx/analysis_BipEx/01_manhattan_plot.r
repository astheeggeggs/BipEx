library(data.table)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

source("r_functions/manhattan_plots.r")

# Common SNPs.
GWAS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.tsv.bgz | gzcat'
gwas <- fread(GWAS_TSV, header=TRUE) %>% 
    mutate(CHR = as.integer(gsub("chr", "", sapply(strsplit(locus, split=':'),`[`,1))),
    	   POS = as.integer(gsub("chr", "", sapply(strsplit(locus, split=':'),`[`,2))))

## Just keep CHR, BP, and P
# Should probably impose MAC filters for this.
gwas$CHR[is.na(gwas$CHR)] <- 23
output <- "../../analysis_plots/manhattan_plots/01_EWAS.pdf"
dir_separate_figures_output <- "../../analysis_plots/manhattan_plots/EWAS/"
dir.create(dir_separate_figures_output, showWarnings=FALSE)
separate_figures_output <- paste0(dir_separate_figures_output, "01_EWAS_")
to_loop_over <- c("BP1", "BP2", "BPNOS", "BPSCZ", "BP", "BP_including_BPSCZ", "BPPSY", "BP_no_PSY")
title_to_loop_over <- c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder NOS",
	"Schizoaffective", "Bipolar Disorder", "Bipolar Disorder including Schizoaffective",
	"Bipolar Disorder with psychosis", "Bipolar Disorder without psychosis")

pdf(output, width=9, height=4)
for(i in 1:length(to_loop_over)) {
	p <- make_manhattan_plot(
		gwas$CHR, gwas$POS,
		unlist(gwas[paste0('logreg.', to_loop_over[i], '.p_value')]),
		gwas$gene_symbol,
		threshold=6,
		title=title_to_loop_over[i],
		buffer=0,
		ggplot_theme=theme_minimal,
		colour_1="#c24100",
		colour_2="#0081c2",
		save_figure=TRUE,
		file=paste0(separate_figures_output, to_loop_over[i]),
		title_size=22)
	print(p)
}
dev.off()

# Common variants
output <- "../../analysis_plots/manhattan_plots/01_common_EWAS.pdf"
dir_separate_figures_output <- "../../analysis_plots/manhattan_plots/EWAS_common/"
dir.create(dir_separate_figures_output, showWarnings=FALSE)
separate_figures_output <- paste0(dir_separate_figures_output, "01_EWAS_common_")
gwas_common <- gwas %>% filter(MAC >= 40)

pdf(output, width=9, height=4)
for(i in 1:length(to_loop_over)) {
	p <- make_manhattan_plot(
		gwas_common$CHR, gwas_common$POS,
		unlist(gwas_common[paste0('logreg.', to_loop_over[i], '.p_value')]),
		gwas_common$gene_symbol,
		threshold=6,
		title=title_to_loop_over[i],
		buffer=0,
		ggplot_theme=theme_minimal,
		colour_1="#c24100",
		colour_2="#0081c2",
		save_figure=TRUE,
		file=paste0(separate_figures_output, to_loop_over[i]),
		title_size=22)
	print(p)
}
dev.off()

## Exclude the non-coding variants
output <- "../../analysis_plots/manhattan_plots/01_EWAS_coding.pdf"
dir_separate_figures_output <- "../../analysis_plots/manhattan_plots/EWAS_coding/"
dir.create(dir_separate_figures_output, showWarnings=FALSE)
separate_figures_output <- paste0(dir_separate_figures_output, "01_EWAS_coding_")
gwas <- gwas %>% filter(consequence_category != "non_coding")
pdf(output, width=9, height=4)
for(i in 1:length(to_loop_over)) {
	p <- make_manhattan_plot(
		gwas$CHR, gwas$POS,
		unlist(gwas[paste0('logreg.', to_loop_over[i], '.p_value')]),
		gwas$gene_symbol,
		threshold=6,
		title=title_to_loop_over[i],
		buffer=0,
		ggplot_theme=theme_minimal,
		colour_1="#c24100",
		colour_2="#0081c2",
		save_figure=TRUE,
		file=paste0(separate_figures_output, to_loop_over[i]),
		title_size=22)
	print(p)
}
dev.off()

# Common variants
output <- "../../analysis_plots/manhattan_plots/01_common_EWAS_coding.pdf"
dir_separate_figures_output <- "../../analysis_plots/manhattan_plots/EWAS_common_coding/"
dir.create(dir_separate_figures_output, showWarnings=FALSE)
separate_figures_output <- paste0(dir_separate_figures_output, "01_EWAS_common_coding_")
gwas_common <- gwas_common %>% filter(consequence_category != "non_coding")

pdf(output, width=9, height=4)
for(i in 1:length(to_loop_over)) {
	p <- make_manhattan_plot(
		gwas_common$CHR, gwas_common$POS,
		unlist(gwas_common[paste0('logreg.', to_loop_over[i], '.p_value')]),
		gwas_common$gene_symbol,
		threshold=6,
		title=title_to_loop_over[i],
		buffer=0,
		ggplot_theme=theme_minimal,
		colour_1="#c24100",
		colour_2="#0081c2",
		save_figure=TRUE,
		file=paste0(separate_figures_output, to_loop_over[i]),
		title_size=22)
	print(p)
}
dev.off()
