library(data.table)

INITIAL_VARIANT_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_metrics.tsv'
INITIAL_VARIANT_LIST <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

dt <- fread(INITIAL_VARIANT_QC_FILE)
n_initial_variants <- nrow(fread(INITIAL_VARIANT_LIST))

summary_df <- data.frame(
	Filter=c(
		"Variants with < 7 alleles",
		"Failing VQSR",
		"In LCRs",
		"Outside padded target interval",
		"Invariant sites after initial variant and genotype filters",
		"Variants after initial filtering"
		),
	Variants=c(
		nrow(dt),
		sum(dt$fail_VQSR),
		sum(dt$in_LCR),
		sum(dt$not_in_padded_target_intervals),
		nrow(dt) - sum(dt$not_in_padded_target_intervals | dt$fail_VQSR | dt$in_LCR) - n_initial_variants,
		n_initial_variants
		)
)

fwrite(summary_df, "../../variants_BipEx/02_summary_variant_table.tsv", row.names=FALSE, sep='\t')
