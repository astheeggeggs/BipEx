library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)

source("r_functions/burden_tests.r")

# Using Gene ID.
## TSV files
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC5_gene_counts_per_sample_ENSG.tsv.bgz | zcat'
GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC5_gene_counts_per_sample_ENSG.tsv.bgz | zcat'
# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz | zcat'
GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz | zcat'

gene_list_information_GTEx <- fread("../../sLDSC_GTEx_geneset_data/GTEx_tissue_information.tsv")$Tissue
for (i in 1:length(gene_list_information_GTEx)) {
	if  (i == 1) {
		gene_list_GTEx <- fread(paste0("../../sLDSC_GTEx_geneset_data/GTEx.", i, ".GeneSet"), header=FALSE) %>% dplyr::rename(gene=V1) %>% mutate(geneset_name=gene_list_information_GTEx[i])
	} else {
		gene_list_GTEx_tmp <- fread(paste0("../../sLDSC_GTEx_geneset_data/GTEx.", i, ".GeneSet"), header=FALSE) %>% dplyr::rename(gene=V1) %>% mutate(geneset_name=gene_list_information_GTEx[i])
		gene_list_GTEx <- rbind(gene_list_GTEx, gene_list_GTEx_tmp)
	}
}

create_genelist_sample_count_tables <- TRUE
consequence_categories <- c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv")
       
# This code reads in count information and gene list information, and combines to create a matrix that is 
# samples x gene_lists, containing the counts in each of the gene lists for each sample.

files_to_read <- c(
	GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV,
	GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV,
	GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV,
	GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV
)

file_out <- c(
	"BP_including_BPSCZ_MAC5_GTEx_gene_set_counts_per_sample.tsv",
	"BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample.tsv",
	"SCZ_MAC5_GTEx_gene_set_counts_per_sample.tsv",
	"SCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample.tsv"
)

if (create_genelist_sample_count_tables)
{
	for (i in 1:length(files_to_read)) {
		dt <- fread(files_to_read[i]) %>% rename(gene_symbol=gene_id)
		dt_to_write <- create_burden_file(dt, gene_list_GTEx, consequence_categories)
		fwrite(dt_to_write, sep='\t', paste0("/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/", file_out[i]))
		# gzip(paste0("gene_lists_output/", file_out[i]), destname=paste0("gene_lists_output/", file_out[i], ".gz"))
	}
}

system("gsutil cp /psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/* gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/")

# Now let's perform the analyses.
# Note, this is running the same code as the permutation tests in pipeline.
# (see dockerfiles/run_observed.r)

gene_set_counts_file_str <- c(
	"/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/BP_including_BPSCZ_MAC5_GTEx_gene_set_counts_per_sample.tsv",
	"/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_output/BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample.tsv"
	)

system("mkdir -p /psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output")


observation_output_coding_burden_files <- c(
	'/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output/BP_including_BPSCZ_MAC5_GTEx_gene_set_counts_per_sample_observed_coding_burden.tsv',
	'/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output/BP_including_BPSCZ_MAC5_GTEx_gnom_non_psych_gene_set_counts_per_sample_observed_coding_burden.tsv'
)

# Move the phenotype information from the cloud so that we can use the code as is
system("gsutil cp gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_burden_BP_including_BPSCZ.tsv.bgz /psych/genetics_data/dpalmer/bipolar_genesets/")
sample_burden_file_str <- '/psych/genetics_data/dpalmer/bipolar_genesets/02_sample_burden_BP_including_BPSCZ.tsv.bgz'

for (i in 1:length(gene_set_counts_file_str)) {
	system(paste("Rscript dockerfiles/run_observed.r", sample_burden_file_str, gene_set_counts_file_str[i], observation_output_coding_burden_files[i], "coding_burden"))
}

# Copied down from the cluster using
# scp dpalmer@login01:/psych/genetics_data/dpalmer/bipolar_genesets/gene_lists_observed_GTEx_output/*observed*tsv Repositories/BipEx/sLDSC_GTEx_geneset_data/gene_lists_observed_GTEx_output/

# Redefine local name for observed burden p-values.
observation_output_overall_burden_files <- dir("../../sLDSC_GTEx_geneset_data/gene_lists_observed_GTEx_output/", full.name=TRUE)

# Create the plots
pdf("../../sLDSC_GTEx_geneset_data/geneset_burden_results/GTEx_burden_plots_BP1.pdf", width=6, height=3)
for (file in observation_output_overall_burden_files) {
	dt_plot <- fread(file)
	dt_plot <- dt_plot %>% filter(!gene_sets %in% c("Uterus", "Vagina", "Ovary", "Testis", "Prostate", "Cervix Endocervix", "Fallopian Tube", "Cervix Ectocervix", "Cells EBV-transformed lymphocytes", "Cells Transformed fibroblasts"))
	dt_plot <- dt_plot %>% mutate(brain = factor(ifelse(grepl("brain", gene_sets, ignore.case=TRUE), 'Brain', 'Other'), levels=c("Brain", "Other")))
	dt_plot <- dt_plot %>% arrange(brain, obs_log_is_BP1_ptv_p.value) %>% mutate(order=seq(1,nrow(dt_plot)))
	p <- ggplot(dt_plot, aes(x=order, y=-log10(obs_log_is_BP1_ptv_p.value), col=brain, fill=brain)) + 
	geom_col() + theme_classic() + ylab(expression(paste(-log[10], '(p-value)'))) + xlab('GTEx tissue (sorted by group, then p-value)') + theme(legend.title = element_blank())
	print(p)
	print(head(dt_plot, 14))
}
dev.off()

pdf("../../sLDSC_GTEx_geneset_data/geneset_burden_results/GTEx_burden_plots_BP2.pdf", width=6, height=3)
for (file in observation_output_overall_burden_files) {
	dt_plot <- fread(file)
	dt_plot <- dt_plot %>% filter(!gene_sets %in% c("Uterus", "Vagina", "Ovary", "Testis", "Prostate", "Cervix Endocervix", "Fallopian Tube", "Cervix Ectocervix", "Cells EBV-transformed lymphocytes", "Cells Transformed fibroblasts"))
	dt_plot <- dt_plot %>% mutate(brain = factor(ifelse(grepl("brain", gene_sets, ignore.case=TRUE), 'Brain', 'Other'), levels=c("Brain", "Other")))
	dt_plot <- dt_plot %>% arrange(brain, obs_log_is_BP2_ptv_p.value) %>% mutate(order=seq(1,nrow(dt_plot)))
	p <- ggplot(dt_plot, aes(x=order, y=-log10(obs_log_is_BP2_ptv_p.value), col=brain, fill=brain)) + 
	geom_col() + theme_classic() + ylab(expression(paste(-log[10], '(p-value)'))) + xlab('GTEx tissue (sorted by group, then p-value)') + theme(legend.title = element_blank())
	print(p)
	print(head(dt_plot, 14))
}
dev.off()

pdf("../../sLDSC_GTEx_geneset_data/geneset_burden_results/GTEx_burden_plots_BP.pdf", width=6, height=3)
for (file in observation_output_overall_burden_files) {
	dt_plot <- fread(file)
	dt_plot <- dt_plot %>% filter(!gene_sets %in% c("Uterus", "Vagina", "Ovary", "Testis", "Prostate", "Cervix Endocervix", "Fallopian Tube", "Cervix Ectocervix", "Cells EBV-transformed lymphocytes", "Cells Transformed fibroblasts"))
	dt_plot <- dt_plot %>% mutate(brain = factor(ifelse(grepl("brain", gene_sets, ignore.case=TRUE), 'Brain', 'Other'), levels=c("Brain", "Other")))
	dt_plot <- dt_plot %>% arrange(brain, obs_log_is_BP_ptv_p.value) %>% mutate(order=seq(1,nrow(dt_plot)))
	p <- ggplot(dt_plot, aes(x=order, y=-log10(obs_log_is_BP_ptv_p.value), col=brain, fill=brain)) + 
	geom_col() + theme_classic() + ylab(expression(paste(-log[10], '(p-value)'))) + xlab('GTEx tissue (sorted by group, then p-value)') + theme(legend.title = element_blank())
	print(p)
	print(head(dt_plot, 14))
}
dev.off()
