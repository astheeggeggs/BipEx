# 00_get_bam_info.r

VCF_SAMPLE_IDS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv'
VCF_META <- '/seq/dax/BiPolar_CasesControls1_Exomes/Exome/v2/BiPolar_CasesControls1_Exomes.calling_metadata.txt'
BAM_METRICS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/bam_metrics.tsv'

# 03_initial_sample_qc_filter.r
QC_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv"
SAMPLE_LIST_INITIAL_QC <- '../../samples_BipEx/03_initial_qc.keep.sample_list'

# 03_initial_sample_qc_plot.r
PLOTS <- '../../QC_plots/sample_plots/'
# Define some thresholds 
T_sample_callRate <- 0.93
# T_pct_contaminination <- 0.004
T_pct_contamination <- 0.02
T_pct_chimeras <- 0.015
T_dpMean <- 30
T_gqMean <- 55

# 05_impute_sex_plot.r
IMPUTESEX_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.tsv"
SEXCHECK_LIST <- '../../samples_BipEx/05_sexcheck.remove.sample_list'
Y_NCALLED_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_ycalled.tsv"
T_impute_sex <- 0.6

# 06_ibd_plot.r
IBD_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.tsv"
IBD_THRESHOLD <- 0.2

# 06_ibd_filtered.r
SAMPLE_LIST_IBD <- '../../samples_BipEx/06_ibd.remove.sample_list'

# 08_ultra_rare_counts_plot.r
URV_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv'

# 09_10_pca_plot.r
PCA_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/09_pca_scores.tsv'
PCA_1KG_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_pca_scores_1kg.tsv'
EUROPEAN_SAMPLES_STRICT <- '../../samples_BipEx/10_european.strict.sample_list'
EUROPEAN_SAMPLES_LOOSE <- '../../samples_BipEx/10_european.loose.sample_list'
EUROPEAN_SAMPLES_EXCLUDING_URV_OUTLIERS <- '../../samples_BipEx/10_european.no_URV_outliers.sample_list'

T_nURVSNP <- 300
T_nURVIndel <- 25

T_European_RF <- 0.95

PCA_EUR_SCORES <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/11_pca_scores.strict_european.tsv"

# 12_pca_us_PCA_plot.r
PCA_SCORES_USA <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_pca_scores_USA_samples.tsv"

EUROPEAN_AND_AJ_SAMPLES <- "../../samples_BipEx/12_european_and_aj.sample_list"
EUROPEANS <- "../../samples_BipEx/12_european.sample_list"
MAINLAND_EUROPEANS <- "../../samples_BipEx/12_mainland_european.sample_list"
SWEDES <- "../../samples_BipEx/12_swedes.sample_list"

# 13_pca_AJ_1kg_plots.r
PCA_EUR_1KG_AJ_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.european_and_aj.1kg.tsv'
PCA_EUR_1KG_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.european.1kg.tsv'

# 14_final_variant_qc_plot.r
VARIANT_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.variants.tsv.bgz | gzcat'
T_variant_call_rate  <- 0.97
T_absdiff <- 0.02
T_pHWE <- 1e-6

# 14_final_variant_qc_filter.r
VARIANT_LIST <- '../../variants_BipEx/14_final_qc.keep.variant_list'

# 15_final_sample_qc_plot.r
SAMPLE_BEFORE_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.after.samples.tsv'

# 15_final_sample_qc_filter.r
FINAL_SAMPLE_LIST <- '../../samples_BipEx/15_final_qc.keep.sample_list'
FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS <- '../../samples_BipEx/15_final_qc_remove_singleton_outliers.keep.sample_list'

# 17_pca_final_plot.r
FINAL_SAMPLE_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz | gzcat'

