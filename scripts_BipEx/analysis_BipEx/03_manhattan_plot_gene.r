library(data.table)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

source("r_functions/manhattan_plots.r")
location_of_p_value_tables <- "../../analysis_plots/gene_counts_qq/plots/"

hgnc <- fread("../browser_BipEx/hgnc.tsv") %>% rename("gene_symbol" = `Approved symbol`)
ensembl <- fread("ensembl_IDs_and_positions.tsv") %>% rename(ensembl_ID = `Gene stable ID`)

titles <- c("BP1", "BP2", "BP")
titles_translation <- c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder")

# Read in the p-value information for each of the results that are so far plotted in QQ plots
for (file in paste0(location_of_p_value_tables, grep("BP_gene_CMH_MAC5_gnom_non_psych_qq.tsv", dir(location_of_p_value_tables), value=TRUE)))
{
	dt <- fread(file)
	dt <- merge(dt, hgnc, by="gene_symbol") %>% rename(ensembl_ID = `Ensembl ID(supplied by Ensembl)`)
	dt <- merge(dt, ensembl, by="ensembl_ID")

	dt <- dt %>% rename(
		CHR = `Chromosome/scaffold name`,
		START = `Gene start (bp)`,
		END = `Gene end (bp)`,
	) %>% mutate(BP = (START + (END-START)/2))

	if (grepl("CMH", file)) {
		dt <- dt %>% select(-c(grep('BPNOS', names(dt), value=TRUE), grep('BPSCZ', names(dt), value=TRUE)))
	}
	pdf(gsub("qq\\.tsv", "manhattan.pdf", file), width=10, height=5)
	for(pval in grep("pval", names(dt), value=TRUE)) {

		title_tmp <- titles_translation[which(titles == gsub(".*is_(.*)\\..*", "\\1", pval))]
		# consequences <- c("non_coding", "synonymous", "damaging_missense", "ptv")
		consequences <- c("ptv")
		# consequences_translation <- c("Non-coding", "Synonymous", "Damaging missense", "PTV")
		consequences_translation <- c("PTV")		

		for(i in 1:length(consequences)) {
			dt_tmp <- dt %>% filter(consequence_category == consequences[i])
			p <- make_manhattan_plot(dt_tmp$CHR, dt_tmp$BP, dt_tmp[[pval]], dt_tmp$gene_symbol,
				log_p_vals=FALSE, significance_T=2.14e-6, threshold=1,
				title=paste0(title_tmp, ": ", consequences_translation[i]),
				colour_1="#c24100",
				colour_2="#0081c2")
			print(p)
		}
	}
	dev.off()
}


# # For the website
# website_plots <- "../../analysis_plots/gene_counts_qq/website_plots/"
# for (file in paste0(location_of_p_value_tables, grep("BP_gene_.*tsv", dir(location_of_p_value_tables), value=TRUE)))
# {
# 	cat(file, '\n')
# 	dt <- fread(file)
# 	dt <- merge(dt, hgnc, by="gene_symbol") %>% rename(ensembl_ID = `Ensembl ID(supplied by Ensembl)`)
# 	dt <- merge(dt, ensembl, by="ensembl_ID")

# 	dt <- dt %>% rename(
# 		CHR = `Chromosome/scaffold name`,
# 		START = `Gene start (bp)`,
# 		END = `Gene end (bp)`,
# 	) %>% mutate(BP = (START + (END-START)/2))

# 	type <- paste0(
# 		ifelse(grepl("CMH", file), 'CMH', 'fisher'),
# 		ifelse(grepl("MAC5", file), '_MAC5', ''),
# 		ifelse(grepl("non_psych", file), '_gnomad_non_psych', '')
# 	)

# 	if (grepl("CMH", file)) {
# 		dt <- dt %>% select(-c(grep('is_BPNOS', names(dt), value=TRUE), grep('is_BPSCZ', names(dt), value=TRUE)))
# 	}

# 	for(pval in grep("pval", names(dt), value=TRUE)) {

# 		title_tmp <- paste0(type, '_', titles[which(titles == gsub(".*is_(.*)\\..*", "\\1", pval))])
# 		cat(title_tmp, '\n')
# 		consequences <- c("non_coding", "synonymous", "damaging_missense", "ptv")
# 		consequences_translation <- c("Non-coding", "Synonymous", "Damaging missense", "PTV")

# 		for(i in 1:length(consequences)) {
# 			dt_tmp <- dt %>% filter(consequence_category == consequences[i])
# 			p <- make_manhattan_plot(
# 				dt_tmp$CHR, dt_tmp$BP,
# 				dt_tmp[[pval]],
# 				dt_tmp$gene_symbol,
# 				log_p_vals=FALSE,
# 				significance_T=2.14e-6,
# 				threshold=3,
# 				title="", #consequences_translation[i],
# 				file=paste0(website_plots, title_tmp, "_", consequences[i], "_manhattan"),
# 				save_figure=TRUE,
# 				title_size=22,
# 				two_tone=TRUE,
# 				colour_1="#c24100",
# 				colour_2="#0081c2",
# 				ggplot_theme=theme_minimal,
# 				buffer=0)
# 		}
# 	}
# }
