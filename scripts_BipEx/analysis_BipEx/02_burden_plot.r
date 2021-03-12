library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(logistf)
library(forestplot)
library(RColorBrewer)
library(wesanderson)

source("../QC_BipEx/r_functions_and_parameters/pretty_plotting.r")
source("r_functions/burden_tests.r")

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
		is_BP_no_PSY = "phenotype_boolean.is_BP_no_PSY"
		)

dt_BP_including_BPSCZ <- dt_BP_including_BPSCZ %>% select(-c(grep("MPC_3", names(dt_BP_including_BPSCZ), value=TRUE)))

dt_BP_including_BPSCZ <- create_forest_names(dt_BP_including_BPSCZ)
dt_BP_including_BPSCZ <- tidy_data(dt_BP_including_BPSCZ)
tests_BP_including_BPSCZ <- get_tests_and_covariates(dt_BP_including_BPSCZ)

load("forest_plot_Rdata_files/BP_including_BPSCZ.Rdata")

# Make sure folders are created.
system("mkdir ../../analysis_plots/forest_plots/small_forest_plot")
system("mkdir ../../analysis_plots/forest_plots/small_forest_plot/BP_including_BPSCZ")
system("mkdir ../../analysis_plots/forest_plots/small_forest_plot/SCZ")

# DEV: For these, include the p-values and sample sizes in the titles.
small_forest_dir <- '../../analysis_plots/forest_plots/small_forest_plot/'
small_forest_BP_dir <- paste0(small_forest_dir, 'BP_including_BPSCZ/')
small_forest_SCZ_dir <- paste0(small_forest_dir, 'SCZ/')

large_forest_dir <- '../../analysis_plots/forest_plots/large_forest_plot/'
large_forest_BP_dir <- paste0(large_forest_dir, 'BP_including_BPSCZ/')
large_forest_SCZ_dir <- paste0(large_forest_dir, 'SCZ/')

clip_OR <- c(0.8, 1.4)
pretitles <- c("Bipolar Disorder 1: ", "Bipolar Disorder 2: ", "Bipolar Disorder: ", "Bipolar Disorder including Schizoaffective: ",
	"Bipolar Disorder NOS: ", "Schizoaffective: ", "Bipolar Disorder with Psychosis: ", "Bipolar Disorder without Psychosis: ")

i <- 1
for (test in names(forest_plots)) {

	test_name <- gsub("is_", "", test)
	create_small_forest_plot(forest_plots_coding_burden[[test]][["lin"]], 'label',
		paste0(small_forest_BP_dir, test_name, "_lin_coding_burden.pdf"), 
		zero=0, xlabel='Excess variants per case', table_cols=c("mean", "p_vals"),
		table_col_names=c("Variants", "Excess", "p-value"), pretitle=pretitles[i])

	create_small_forest_plot(forest_plots_coding_burden[[test]][['logit']], 'label',
		paste0(small_forest_BP_dir, test_name, "_logit_coding_burden.pdf"), 
		scale_boxes=0.5, xlabel='Odds ratio', table_cols=c("mean", "p_vals"),
		table_col_names=c("Variants", "OR", "p-value"), clip=clip_OR,
		pretitle=pretitles[i])
	i <- i+1
}

load("forest_plot_Rdata_files/SCZ.Rdata")

# Controlling for coding burden
create_small_forest_plot(forest_SCZ_coding_burden$logit, 'label',
	paste0(small_forest_SCZ_dir, 'SCZ_logit_coding_burden.pdf'),
	xlabel='Odds ratio', table_cols=c("mean", "p_vals"),
	table_col_names=c("Variants", "OR", "p-value"), clip=clip_OR,
	pretitle="Schizophrenia: ")
create_small_forest_plot(forest_SCZ_coding_burden$lin, 'label',
	paste0(small_forest_SCZ_dir, 'SCZ_linear_coding_burden.pdf'),
	zero=0, xlabel='Excess variants per case', table_cols=c("mean", "p_vals"),
	table_col_names=c("Variants", "Excess", "p-value"), pretitle="Schizophrenia: ")

load("forest_plot_Rdata_files/BP_including_BPSCZ_locations.Rdata")
system("mkdir ../../analysis_plots/forest_plots/large_forest_plot/BP_including_BPSCZ")

i <- 1
for (test in names(forest_plots_locations))
{
	cat(test, "...")
	test_name <- gsub("is_", "", test)
	for (model in c("lin", "logit"))
	{	
		cat(model, "...")

		full_dt_tmp_coding_burden <- cbind(
			data.frame(variants=rep(tests_BP_including_BPSCZ$tests[1],
				nrow(forest_plots_locations_coding_burden[[test]][[model]][[1]]))),
			forest_plots_locations_coding_burden[[test]][[model]][[1]]
		)

		for(j in 2:length(tests_BP_including_BPSCZ$tests))
		{
			full_dt_tmp_coding_burden <- rbind(
				full_dt_tmp_coding_burden, cbind(
					data.frame(variants=rep(tests_BP_including_BPSCZ$tests[j],
						nrow(forest_plots_locations_coding_burden[[test]][[model]][[j]]))),
					forest_plots_locations_coding_burden[[test]][[model]][[j]]
				)
			)
		}

		if(model == "logit") {
			create_large_forest_plot(full_dt_tmp_coding_burden,
				paste0(large_forest_BP_dir, test_name, "_split_by_location_", model, "_coding_burden.pdf"),
				tests_BP_including_BPSCZ$tests, xlabel='Odds ratio', clip=clip_OR, pretitle=pretitles[i])
		} else {
			create_large_forest_plot(full_dt_tmp_coding_burden,
				paste0(large_forest_BP_dir, test_name, "_split_by_location_", model, "_coding_burden.pdf"),
				tests_BP_including_BPSCZ$tests, zero=0, xlabel='Excess variants per case', pretitle=pretitles[i])
		}
	}
	cat('\n')
	i <- i+1
}

