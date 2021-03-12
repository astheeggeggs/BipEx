library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(logistf)
library(forestplot)
library(RColorBrewer)
library(wesanderson)
library(xtable)

source("r_functions/burden_tests.r")
rerun <- FALSE

PLOT <- "../../analysis_plots/"
SAMPLE_BURDEN_FILE_BP_including_BPSCZ = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_burden_BP_including_BPSCZ.tsv.bgz | gzcat'
SAMPLE_BURDEN_FILE_SCZ = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_burden_SCZ.tsv.bgz | gzcat'

dt_BP_including_BPSCZ <- fread(SAMPLE_BURDEN_FILE_BP_including_BPSCZ) %>%
	rename(
		is_BP1 = "phenotype_boolean.is_BP1",
		is_BP2 = "phenotype_boolean.is_BP2",
		is_BP = "phenotype_boolean.is_BP",
		is_BP_including_BPSCZ = "phenotype_boolean.is_BP_including_BPSCZ",
		is_BPNOS = "phenotype_boolean.is_BPNOS",
		is_BPSCZ = "phenotype_boolean.is_BPSCZ",
		is_BPPSY = "phenotype_boolean.is_BPPSY",
		is_BP_no_PSY = "phenotype_boolean.is_BP_no_PSY",
		PSYCHOSIS = "phenotype.PSYCHOSIS",
		AGE_FI_24 = "phenotype.AGE_FI_24",
		AGE_FI_40 = "phenotype.AGE_FI_40",
		AGE_FS_24 = "phenotype.AGE_FS_24",
		AGE_FS_40 = "phenotype.AGE_FS_40",
		AGE_D_24 = "phenotype.AGE_D_24",
		AGE_D_40 = "phenotype.AGE_D_40"
		)

# Do not include controls from mismatched ancestry.
dt_BP_including_BPSCZ <- dt_BP_including_BPSCZ %>%
	filter(phenotype.LOCATION %in% c("Boston, USA", "Cardiff, UK", "London, UK", "Stockholm, SWE", "Wurzburg, GER"))

# Exclude MPC > 3. Not sufficient samples in the smaller cohorts.
dt_BP_including_BPSCZ <- dt_BP_including_BPSCZ %>% select(-c(grep("MPC_3", names(dt_BP_including_BPSCZ), value=TRUE)))

cat("\nCreating groupings of locations for forest plots\n")
dt_BP_including_BPSCZ <- create_forest_names(dt_BP_including_BPSCZ)

Locations <- names(table(dt_BP_including_BPSCZ$Forest_location))
count_tests <- names(dt_BP_including_BPSCZ)[grep("burden", names(dt_BP_including_BPSCZ))]
# Simplify things by just looking at the overall burden (counting both SNPs and indels).
where <- c(grep("indel", count_tests), grep("SNP", count_tests))
count_tests <- count_tests[-where]

phenotypes_to_test <- c("is_BP1", "is_BP2", "is_BP")

## Regressions, and meta-analysis thereof.
tests_BP_including_BPSCZ <- get_tests_and_covariates(dt_BP_including_BPSCZ)
dt_BP_including_BPSCZ <- tidy_data(dt_BP_including_BPSCZ)

if (!file.exists("forest_plot_Rdata_files/BP_including_BPSCZ_psychosis.Rdata") | (rerun == TRUE)) {
	
	forest_plots_coding_burden <- list()

	for(test in phenotypes_to_test)
	{
		# Including total coding as a covariate.
		cat("Testing boolean phenotype", test, "with n coding burden as a covariate...\n")
		forest_plots_coding_burden[[test]] <- run_collection_burden_regression(
			tests_BP_including_BPSCZ$tests,
			dt_BP_including_BPSCZ %>% filter(dt_BP_including_BPSCZ[[test]] == TRUE) %>% filter(!is.na(PSYCHOSIS)), 'PSYCHOSIS',
			run_linear=TRUE, run_poisson=FALSE, run_logistic=FALSE, count_tests_covariates=tests_BP_including_BPSCZ$covariates_coding
		)
	}

	# Save the results to an Rdata file for plotting.
	save(forest_plots_coding_burden, file = "forest_plot_Rdata_files/BP_including_BPSCZ_psychosis.Rdata")
}

cat("Now testing location by location and performing meta-analysis...\n\n")
