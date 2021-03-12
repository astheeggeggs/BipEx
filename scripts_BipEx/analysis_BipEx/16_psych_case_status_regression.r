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

dt_BP_including_BPSCZ <- tidy_data(dt_BP_including_BPSCZ)
dt <- dt_BP_including_BPSCZ %>% filter(!is.na(phenotype.PSYCHOSIS) | !is_BP) %>% mutate(PSYCHOSIS = ifelse(is.na(phenotype.PSYCHOSIS), FALSE, phenotype.PSYCHOSIS))
test_burden(dt, "burden_gnom_non_psych_pli_09.n_URV_PTV", "is_BP", covariates=c("PSYCHOSIS", "burden_gnom_non_psych_pli_09.n_coding_URV"))
