library(data.table)
library(dplyr)

# In this script we pace through the bins - there are 6 bins that we can construct from the age at first impairment data. 

# We will find the case-case comparison with the lowest p-value, and then perform shuffles to evaluate an empirical p-value.

# Read in the data and restrict to the samples for which we have age of onset information.
SAMPLE_BURDEN_FILE_BP_including_BPSCZ = 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_burden_BP_including_BPSCZ.tsv.bgz | gzcat'

dt <- fread(SAMPLE_BURDEN_FILE_BP_including_BPSCZ)
dt <- dt %>% filter(!is.na(phenotype.AGE_FI_24) & !is.na(phenotype.AGE_FI_40))

# Define 6 bins

dt <- dt %>% mutate(AGE_FI_5_bins = 
	case_when(
		(phenotype.AGE_FI_24 == 0) ~ 0,
		((phenotype.AGE_FI_24 == 1) & (phenotype.AGE_FI_40 == 0)) ~ 1,
		((phenotype.AGE_FI_24 == 1) & (phenotype.AGE_FI_40 == 1)) ~ 2,
		((phenotype.AGE_FI_24 == 2) & (phenotype.AGE_FI_40 == 1)) ~ 3,
		((phenotype.AGE_FI_24 == 2) & (phenotype.AGE_FI_40 == 2)) ~ 4,
		TRUE ~ NA_real_
		)
	)		

# Now, check the singleton ptv counts across these and run Kolmogorov-Smirnov.

dt <- dt %>% mutate(
	ptv_carrier = ifelse(burden.n_URV_PTV > 0, TRUE, FALSE),
	ptv_gnom_non_psych_carrier = ifelse(burden_gnom_non_psych.n_URV_PTV > 0, TRUE, FALSE),
	ptv_gnom_non_psych_pli_09_carrier = ifelse(burden_gnom_non_psych_pli_09.n_URV_PTV > 0, TRUE, FALSE),
	ptv_gnom_non_psych_pli_0995_carrier = ifelse(burden_gnom_non_psych_pli_0995.n_URV_PTV > 0, TRUE, FALSE)
	)

# For each, step through the age categories and find the most extreme.

younger_group <- c(0, 1, 2, 3)
older_group <- c(4, 3, 2, 1)

data.table(younger_group=NA, older_group=NA)

counts_to_check <- c(
	"ptv_carrier",
	"ptv_gnom_non_psych_carrier",
	"ptv_gnom_non_psych_pli_09_carrier",
	"ptv_gnom_non_psych_pli_0995_carrier"
	)

burden_to_check <- c(
	"burden.n_URV_PTV",
	"burden_gnom_non_psych.n_URV_PTV",
	"burden_gnom_non_psych_pli_09.n_URV_PTV",
	"burden_gnom_non_psych_pli_0995.n_URV_PTV"
	)

min_KS_p <- Inf
min_Fisher_p <- Inf

for(k in 1:length(counts_to_check))
{
	count <- counts_to_check[k]
	burden <- burden_to_check[k]

	for(i in 1:length(younger_group))
	{
		younger_group_tmp <- younger_group[1:i]
		for(j in 1:(length(older_group) - length(younger_group_tmp) + 1))
		{
			older_group_tmp <- older_group[1:j]
			cat("current comparison:\n")
			cat(paste0(count, "\n"))
			cat(paste(younger_group_tmp, collapse=","), "against", paste(sort(older_group_tmp), collapse=","), "\n\n")
			
			dt <- dt %>% mutate(
				current_test = case_when(
					AGE_FI_5_bins %in% younger_group_tmp ~ "younger",
					AGE_FI_5_bins %in% older_group_tmp ~ "older",
					TRUE ~ NA_character_))

			younger_group_carrier_status <- as.integer(unlist(dt %>% filter(current_test == "younger") %>% select(!!count)))
			older_group_carrier_status <- as.integer(unlist(dt %>% filter(current_test == "older") %>% select(!!count)))
			

			younger_group_burden <- as.integer(unlist(dt %>% filter(current_test == "younger") %>% select(!!burden)))
			older_group_burden <- as.integer(unlist(dt %>% filter(current_test == "older") %>% select(!!burden)))

			ks_current <- ks.test(
				younger_group_carrier_status,
				older_group_carrier_status
				)
			cat(paste0("KS p-value: ", (ks_current$p.value)), "\n")

			ks_current <- ks.test(
				younger_group_burden, 
				older_group_burden
				)
			cat(paste0("KS p-value: ", (ks_current$p.value)), "\n")

			contingency <- matrix(
					c(sum(younger_group_carrier_status == 1), sum(younger_group_carrier_status == 0),
					  sum(older_group_carrier_status == 1), sum(older_group_carrier_status == 0)),
					nrow = 2)

			print(contingency)
			fisher_current <- fisher.test(contingency)
			cat(paste0("Fisher p-value: ", fisher_current$p.value), "\n")

			min_KS_p <- min(min_KS_p, ks_current$p.value)
			min_Fisher_p <- min(min_Fisher_p, fisher_current$p.value)

		}
	}
}
