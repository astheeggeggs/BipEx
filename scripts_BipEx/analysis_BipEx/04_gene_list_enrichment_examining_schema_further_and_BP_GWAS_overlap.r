library(data.table)
library(dplyr)
library(R.utils)

source("r_functions/burden_tests.r")

# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz | zcat'

# Looking further at SCHEMA, and SCZ GWAS overlap.
# SCZ finemapping is in SCHEMA Table S11.

gene_lists_schema <- fread("../../SCHEMA_data/meta_results_2020_10_05_16_35_58.txt") %>% head(200) %>% 
mutate(gene_set = c(
	rep("schema_1_50", 50),
	rep("schema_51_100", 50),
	rep("schema_101_150", 50), 
	rep("schema_151_200", 50)
	)
) %>% select(Gene, gene_set) %>% rename(gene=Gene, geneset_name=gene_set)

# Next add the SCZ fine-mapped genes.
gene_lists_SCZ <- fread("../../SCHEMA_data/Table_S11.txt") %>% 
	filter(Name == "PGC3 FINEMAP-prioritized genes") %>% 
	select(c("Gene ID", "Name")) %>% rename(gene=`Gene ID`, geneset_name=Name)

gene_lists <- rbind(gene_lists_schema, gene_lists_SCZ)

create_genelist_sample_count_tables <- TRUE
consequence_categories <- rev(c("synonymous", "other_missense", "damaging_missense", "ptv"))

# This code reads in count information and gene list information, and combines to create a matrix that is 
# samples x gene_lists, containing the counts in each of the gene lists for each sample.

file_to_read <- GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV
file_out_prefix <- "BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_SCZ_overlap"

if (create_genelist_sample_count_tables) {
	dt <- fread(file_to_read) %>% rename(gene_symbol=gene_id)
	create_burden_file(dt, gene_lists, consequence_categories, file_out_prefix=file_out_prefix)
	create_exome_burden_file(
		dt, consequence_categories=c("synonymous", "other_missense", "damaging_missense", "ptv"),
		file_out_prefix=file_out_prefix)
	create_exome_burden_file(
		dt, consequence_categories=c("synonymous"),
		file_out_prefix=paste0(file_out_prefix, "_synonymous"))
}

# Merge the burden files, so that they are in the same format as before, 
# zip them and push them to the cloud and then run the regressions in hail using the 
# `run_regressions_obs' function (in 04_gene_list_enrichment_tests.py) in hail.
init <- TRUE
i <- 1
dt_list <- list()
for (consequence in consequence_categories) {
	cat(consequence, "\n")
	if (init) {
		cat("Reading in...\n")
		dt <- fread(paste0(file_out_prefix[i], '_', consequence, '_gene_lists.tsv'))
		dt_list[[i]] <- dt
		cat("Transposing...\n")
		dt <- transpose(dt, make.names="sample", keep.names="gene_set")
		dt[, consequence_category:=consequence]
		init <- FALSE
	} else {
		cat("Reading in...\n")
		dt_tmp <- fread(paste0(file_out_prefix[i], '_', consequence, '_gene_lists.tsv'))
		dt_list[[i]] <- dt_tmp
		cat("Transposing...\n")
		dt_tmp <- transpose(dt_tmp, make.names="sample", keep.names="gene_set")
		dt_tmp[, consequence_category:=consequence]
		dt <- rbind(dt, dt_tmp)
	}
}

colnames <- names(dt)
new_colnames <- c(colnames[1], colnames[length(colnames)], colnames[2:(length(colnames)-1)])
setcolorder(dt, new_colnames)

fwrite(dt, file=paste0(file_out_prefix[i], '_gene_lists.tsv'), sep='\t')
system(paste0("gsutil cp ", file_out_prefix[i], "_gene_lists.tsv ", "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/"))
system(paste0("gsutil cp ", file_out_prefix[i], "*_burden.tsv ", "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/"))

