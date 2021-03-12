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
SAMPLE_BURDEN_FILE_BP_including_BPSCZ = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_MAC5_burden_BP_including_BPSCZ_no_top_20.tsv.bgz | gzcat'

dt_BP_including_BPSCZ <- fread(SAMPLE_BURDEN_FILE_BP_including_BPSCZ) %>%
	rename(
		is_BP1 = "phenotype_boolean.is_BP1",
		is_BP2 = "phenotype_boolean.is_BP2",
		is_BP = "phenotype_boolean.is_BP",
		is_BP_including_BPSCZ = "phenotype_boolean.is_BP_including_BPSCZ",
		is_BPNOS = "phenotype_boolean.is_BPNOS",
		is_BPSCZ = "phenotype_boolean.is_BPSCZ",
		is_BPPSY = "phenotype_boolean.is_BPPSY",
		is_BP_no_PSY = "phenotype_boolean.is_BP_no_PSY"
		)

# Exclude MPC > 3. Not sufficient samples in the smaller cohorts.
dt_BP_including_BPSCZ <- dt_BP_including_BPSCZ %>% select(-c(grep("MPC_3", names(dt_BP_including_BPSCZ), value=TRUE)))

cat("\nCreating groupings of locations for forest plots\n")
dt_BP_including_BPSCZ <- create_forest_names(dt_BP_including_BPSCZ)

Locations <- names(table(dt_BP_including_BPSCZ$Forest_location))
count_tests <- names(dt_BP_including_BPSCZ)[grep("burden", names(dt_BP_including_BPSCZ))]
# Simplify things by just looking at the overall burden (counting both SNPs and indels).
where <- c(grep("indel", count_tests), grep("SNP", count_tests))
count_tests <- count_tests[-where]

## CMH tests.
phenotypes_to_test <- c("is_BP1", "is_BP2", "is_BP", "is_BPPSY", "is_BP_no_PSY")

if (rerun == TRUE) {
	
	cat("\nBipolar including schizoaffective disorder.\n")
	for(test in phenotypes_to_test) {
		dt_CMH <- run_collection_CMH(dt_BP_including_BPSCZ, count_tests, test)
		xtable(dt_CMH, math.style.exponents = TRUE, auto=TRUE)

	}
}

## Regressions, and meta-analysis thereof.
tests_BP_including_BPSCZ <- get_tests_and_covariates(dt_BP_including_BPSCZ)
dt_BP_including_BPSCZ <- tidy_data(dt_BP_including_BPSCZ)

if (!file.exists("forest_plot_Rdata_files/BP_including_BPSCZ_no_top_20_feb21.Rdata") | (rerun == TRUE)) {
	forest_plots <- list()
	forest_plots_coding_burden <- list()
	forest_plots_no_overall_burden <- list()

	for(test in phenotypes_to_test)
	{
		cat("Testing boolean phenotype", test, "with n coding burden as a covariate...\n")
		forest_plots_coding_burden[[test]] <- run_collection_burden_regression(
			tests_BP_including_BPSCZ$tests,
			dt_BP_including_BPSCZ, test,
			tests_BP_including_BPSCZ$covariates_coding
		)
	}

	# Save the results to an Rdata file for plotting.
	save(forest_plots_coding_burden, file = "forest_plot_Rdata_files/BP_including_BPSCZ_no_top_20_feb21.Rdata")
} else {
	load("forest_plot_Rdata_files/BP_including_BPSCZ_no_top_20_feb21.Rdata")
}

cat("Now testing location by location and performing meta-analysis...\n\n")

forest_plots_locations_coding_burden <- list()

if (!file.exists("forest_plot_Rdata_files/BP_including_BPSCZ_locations_no_top_20_feb21.Rdata") | (rerun == TRUE))
{
	for(test in phenotypes_to_test)
	{
		# Including the covariates (overall burden).
		dt_to_test <- dt_BP_including_BPSCZ

		if (test == "is_BPSCZ") {
			dt_to_test <- dt_to_test %>%
				filter(Forest_location %in% c("NED", "UK/Ireland", "USA"))
		}

		if(test %in% c("is_BPNOS", "is_BPPSY", "is_BP_no_PSY")) {
			dt_to_test <- dt_to_test %>% 
				filter(Forest_location %in% c("GER", "SWE, Stockholm", "UK/Ireland", "USA"))
		}

		cat("Testing boolean phenotype", test, "with n coding burden as a covariate...\n")
		forest_plots_locations_coding_burden[[test]] <- run_locations_regression(
			dt_to_test, tests_BP_including_BPSCZ$tests,
			test, tests_BP_including_BPSCZ$covariates_coding
			)
	}

	save(forest_plots_locations_coding_burden,
		file="forest_plot_Rdata_files/BP_including_BPSCZ_locations_no_top_20_feb21.Rdata")
}

