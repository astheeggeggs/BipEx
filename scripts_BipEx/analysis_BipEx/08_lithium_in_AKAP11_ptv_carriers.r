library(data.table)
library(dplyr)

## TSV files
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_SINGLETON_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_singleton_gene_counts_per_sample.tsv.bgz | gzcat'
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC5_gene_counts_per_sample.tsv.bgz | gzcat'

# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_SINGLETON_GNOM_NON_PSYCH_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_singleton_gnom_non_psych_gene_counts_per_sample.tsv.bgz | gzcat'
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz | gzcat'

files <- c(
	# GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_TSV,
	# GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_SINGLETON_COUNTS_TSV,
	GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV #,
	# GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_SINGLETON_GNOM_NON_PSYCH_COUNTS_TSV
	)

for(file in files) {

	dt <- fread(file)

	dt_AKAP <- dt %>% filter(gene_symbol == "AKAP11") %>% filter(consequence_category == "ptv")

	dt_AKAP <- t(dt_AKAP %>% select(-gene_symbol, -consequence_category))
	dt_AKAP <- data.table(SAMPLE_ALIAS=rownames(dt_AKAP), AKAP_ptv=dt_AKAP[,1])
	dt_AKAP <- dt_AKAP %>% filter(AKAP_ptv != 0)

	# Read in the full phenotype information, and merge with all the available data containing information about Lithium
	dt_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv")

	dt_AKAP <- merge(dt_pheno, dt_AKAP, by="SAMPLE_ALIAS")

	dt_pheno <- fread("../../phenotype_data/ICCBD_data/phenotype_merge_v4.txt") %>% rename(PARTICIPANT_ID=participant_id) %>% select(PARTICIPANT_ID, LI, migraine)
	dt_AKAP <- merge(dt_AKAP, dt_pheno, by="PARTICIPANT_ID", all.x=TRUE)

	print(dt_AKAP)

	print(ks.test((dt_AKAP %>% filter(LI %in% c(0,1,2)))$LI, (dt_pheno %>% filter(LI %in% c(0,1,2)))$LI))
	rm(dt)

}

# Add in the two new good responders that we have from the cardiff data.
print(ks.test(c((dt_AKAP %>% filter(LI %in% c(0,1,2)))$LI, 2), (dt_pheno %>% filter(LI %in% c(0,1,2)))$LI))