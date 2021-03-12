library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)

source("r_functions/burden_tests.r")

# Using Gene ID.
# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz | zcat'

gene_list_information_GTEx <- fread("../../sLDSC_GTEx_geneset_data/GTEx_tissue_information.tsv")$Tissue

for (i in 1:length(gene_list_information_GTEx)) {
	cat(gene_list_information_GTEx[i],"\n")
	if  (i == 1) {
		gene_list_GTEx <- fread(paste0("../../sLDSC_GTEx_geneset_data/GTEx.", i, ".GeneSet"), header=FALSE) %>% dplyr::rename(gene=V1) %>% mutate(geneset_name=gene_list_information_GTEx[i])
	} else {
		gene_list_GTEx_tmp <- fread(paste0("../../sLDSC_GTEx_geneset_data/GTEx.", i, ".GeneSet"), header=FALSE) %>% dplyr::rename(gene=V1) %>% mutate(geneset_name=gene_list_information_GTEx[i])
		gene_list_GTEx <- rbind(gene_list_GTEx, gene_list_GTEx_tmp)
	}
}

gene_list_GTEx <- gene_list_GTEx %>% filter(!geneset_name %in% c("Uterus", "Vagina", "Ovary", "Testis", "Prostate", "Cervix Endocervix", "Fallopian Tube", "Cervix Ectocervix", "Cells EBV-transformed lymphocytes", "Cells Transformed fibroblasts"))

# Gene-sets
brain_genesets <- grep("Brain", unique(gene_list_GTEx$geneset_name), value=TRUE)
gene_list_GTEx <- gene_list_GTEx %>% mutate(geneset_name_brain_non_brain = ifelse(geneset_name %in% brain_genesets, "Brain", "Non-brain"))

# Combine together the brain and non-brain genesets into two new sets, and discard the rest.
gene_list_GTEx <- rbind(
	data.table(gene = unique((gene_list_GTEx %>% filter(geneset_name_brain_non_brain == "Brain"))$gene), geneset_name = "Brain"),
	data.table(gene = unique((gene_list_GTEx %>% filter(geneset_name_brain_non_brain == "Non-brain"))$gene), geneset_name = "Non-brain")
	)

create_genelist_sample_count_tables <- TRUE
consequence_categories <- c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv")
       
# This code reads in count information and gene list information, and combines to create a matrix that is 
# samples x gene_lists, containing the counts in each of the gene lists for each sample.

files_to_read <- c(GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV)

file_out <- "BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample_overall_brain_non_brain.tsv"

if (create_genelist_sample_count_tables)
{
	dt <- fread(files_to_read) %>% rename(gene_symbol=gene_id)
	dt_to_write <- create_burden_file(dt, gene_list_GTEx, consequence_categories, )
	fwrite(dt_to_write, sep='\t', paste0("/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/", file_out))
}

system("gsutil cp /psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/* gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/")

# Now let's perform the analyses.
# Note, this is running the same code as the permutation tests in pipeline.
# (see dockerfiles/run_observed.r)

gene_set_counts_file_str <- "/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample_overall_brain_non_brain.tsv"

system("mkdir -p /psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output")
observation_output_coding_burden_file <- "/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output/BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample_observed_brain_non_brain_coding_burden.tsv"

# Move the phenotype information from the cloud so that we can use the code as is
system("gsutil cp gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_burden_BP_including_BPSCZ.tsv.bgz /psych/genetics_data/dpalmer/bipolar_genesets/")
sample_burden_file_str <- '/psych/genetics_data/dpalmer/bipolar_genesets/02_sample_burden_BP_including_BPSCZ.tsv.bgz'

system(paste("Rscript dockerfiles/run_observed.r", sample_burden_file_str, gene_set_counts_file_str, observation_output_coding_burden_file, "coding_burden"))

# Print out the p-value information for the PTVs.
dt <- fread(observation_output_coding_burden_file)
dt_BP <- dt %>% select(c("gene_sets", grep("obs_log_is_BP_ptv", names(dt), value=TRUE)))
exp(dt_BP$obs_log_is_BP_ptv_coef)
dt_BP <- dt %>% select(c("gene_sets", grep("obs_log_is_BPPSY_ptv", names(dt), value=TRUE)))
exp(dt_BP$obs_log_is_BPPSY_ptv_coef)
