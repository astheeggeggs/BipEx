library(data.table)
library(dplyr)
library(R.utils)

source("r_functions/burden_tests.r")

## TSV files
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC5_gene_counts_per_sample.tsv.bgz | zcat'
GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC5_gene_counts_per_sample.tsv.bgz | zcat'
# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz | zcat'
GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz | zcat'

gene_list_howrigan <- fread("gene_lists_input/lists_howrigan.tsv", sep='\t', header=TRUE)
gene_list_singh <- fread("gene_lists_input/lists_singh.tsv", sep='\t', header=TRUE)
gene_list_eli <- fread("gene_lists_input/list_eli_stahl.tsv", sep='\t', header=TRUE)
gene_list_schema <- fread("gene_lists_input/list_schema.tsv", sep='\t', header=TRUE)

names(gene_list_schema) <- names(gene_list_eli) <- names(gene_list_singh) <- names(gene_list_howrigan)

gene_list <- rbind(gene_list_schema, gene_list_singh, gene_list_eli)
create_genelist_sample_count_tables <- TRUE
consequence_categories <- rev(c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv"))

# This code reads in count information and gene list information, and combines to create a matrix that is 
# samples x gene_lists, containing the counts in each of the gene lists for each sample.
files_to_read <- c(
	GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV,
	GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV
)

file_out_prefix <- c(
	"BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample",
	"SCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample"
)

if (create_genelist_sample_count_tables) {
	for (i in 1:length(files_to_read)) {
		dt <- fread(files_to_read[i])
		create_exome_burden_file(
			dt, consequence_categories=c("synonymous", "other_missense", "damaging_missense", "ptv"),
			file_out_prefix=file_out_prefix[i])
	}
}

# Merge the burden files, so that they are in the same format as before, 
# zip them and push them to the cloud and then run the regressions in hail using the 
# `run_regressions_obs' function (in 04_gene_list_enrichement.py) in hail.
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
dt_sum <- rbindlist(dt_list)[, lapply(.SD, sum, na.rm=TRUE), by=sample]
# This is actually tricky and will lead to the same issue as before.
# As an alternative, let's control for the total burden of the class of variation
# genome-wide, rather than just in the gene-set.
colnames(dt_sum)[-1] <- paste0(colnames(dt_sum)[-1], '_sum')
dt_sum <- transpose(dt_sum, make.names="sample", keep.names="gene_set")

colnames <- names(dt)
new_colnames <- c(colnames[1], colnames[length(colnames)], colnames[2:(length(colnames)-1)])
setcolorder(dt, new_colnames)

fwrite(dt, file=paste0(file_out_prefix[i], '_gene_lists.tsv'), sep='\t')
system(paste0("gsutil cp ", file_out_prefix[i], "_gene_lists.tsv ", "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/"))
system(paste0("gsutil cp ", file_out_prefix[i], "*_burden.tsv ", "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/"))
