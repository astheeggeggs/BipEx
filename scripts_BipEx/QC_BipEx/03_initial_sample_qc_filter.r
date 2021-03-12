library(dplyr)
library(data.table)
source("r_functions_and_parameters/r_options_BipEx.r")

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_sample_qc_plot.r")

# Remove the low coverage samples as defined by Excel spreadsheet e-mailed from Laura Gauthier
low_coverage <- fread("Dalio_Low_Coverage_Samples.txt")$SAMPLE

QC_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv"
df <- fread(QC_FILE, sep='\t', header=TRUE, data.table=FALSE)
names(df) <- gsub("qc\\.", "", names(df))

# Also want to count the number of samples not present in the vcf file.
# Note that throughout here, we are merging in BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv
# rather than relying on the annotated phenotype information, due to changes that occurred.

# df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_final.tsv") %>% rename(s=SAMPLE_ALIAS)

# Updated to correct psychosis and include AAO variables (April 2020)
df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv") %>% rename(s=SAMPLE_ALIAS)

if(any(duplicated(df_pheno$s))) {
    cat('remove duplicated names...\n')
    df_pheno <- df_pheno[-which(duplicated(df_pheno$s)),]
}

df_merged <- merge(df_pheno, df,  by='s')

df_initial_summary_count <- data.table(
	"Filter" = c(
		"Initial samples in vcf", 
		"Unable to obtain both phenotype and sequence information",
		"Unknown phenotype",
		"Low coverage or high contamination"),
	"Samples" = c(nrow(df),
				nrow(df) - nrow(df_merged),
				nrow(filter(df_merged, PHENOTYPE_COARSE == 'Unknown')),
				nrow(df_merged[which(df_merged$s %in% low_coverage),])),
	"Bipolar cases" = c(nrow(df_merged %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				NA,
				NA,
				nrow(df_merged[which((df_merged$s %in% low_coverage) & df_merged$PHENOTYPE_COARSE == "Bipolar Disorder"),])),
	"Controls" = c(nrow(df_merged %>% filter(PHENOTYPE_COARSE == "Control")),
				NA,
				NA,
				nrow(df_merged[which((df_merged$s %in% low_coverage) & df_merged$PHENOTYPE_COARSE == "Control"),])))

df_merged <- df_merged[-which(df_merged$s %in% low_coverage),] %>% filter(PHENOTYPE_COARSE != 'Unknown')
df_initial_summary_count <- rbindlist(list(df_initial_summary_count, list("Samples after initial filter", nrow(df_merged),
	nrow(df_merged %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
	nrow(df_merged %>% filter(PHENOTYPE_COARSE == "Control")))), use.names=FALSE)
fwrite(df_initial_summary_count, file='../../samples_BipEx/03_initial_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

df <- df_merged
names(df) <- gsub("qc_padded_ice\\.", "", names(df))

df_out <- filter(df, call_rate > T_sample_callRate) %>%
	filter(PCT_CONTAMINATION < T_pct_contamination) %>%
	filter(PCT_CHIMERAS < T_pct_chimeras) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
df_summary_count <- data.table(
	"Filter" = c("Samples after initial filter",
			   paste0("Sample call rate < ", T_sample_callRate),
			   paste0("% FREEMIX contamination > ", T_pct_contamination),
			   paste0("% chimeric reads > ", T_pct_chimeras),
			   paste0("Mean DP < ", T_dpMean),
			   paste0("Mean GQ < ", T_gqMean),
			   "Samples after sample QC filters"),
	"Samples" = c(nrow(df),
			    nrow(filter(df, call_rate <= T_sample_callRate)),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination)),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras)),
				nrow(filter(df, dp_stats.mean <= T_dpMean)),
				nrow(filter(df, gq_stats.mean <= T_gqMean)),
				nrow(df_out)),
	"Bipolar Cases" = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				length(which(df_out$s %in% (df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder"))$s))),
	"Controls" = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Control")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				length(which(df_out$s %in% (df %>% filter(PHENOTYPE_COARSE == "Control"))$s))))

fwrite(df_summary_count, file='../../samples_BipEx/03_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
