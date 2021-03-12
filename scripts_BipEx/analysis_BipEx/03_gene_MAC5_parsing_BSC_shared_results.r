library(data.table)
library(dplyr)

# Two stages, we first wish to check that the overall burden of singletons and MAC5 <= 5 is similar to that
# seen in the data analysed so far.

# For these variants, we need to determine which are not in gnomAD, and also intersect them with our variant set.
# We first need to lift them over to build 38 to do this.

# First let's combine the data to more easily import it into hail.
consequence_categories <- c("ptv", "damaging", "other_missense", "synonymous")

init <- TRUE
for (k in 1:length(consequence_categories)) {
	bsc_file <- paste0("../../BSC_data/MAC5/BSC_MAC5.with_ref.", consequence_categories[k])
	if (init) {
		dt <- fread(bsc_file) %>% 
		mutate(
			consequence_category = ifelse(
				consequence_categories[k]=="damaging", "damaging_missense", consequence_categories[k])
		)
		init <- FALSE
	} else {
		dt <- rbind(dt, fread(bsc_file) %>% mutate(consequence_category = consequence_categories[k]))
	}
}

dt <- dt %>% rename(gene_symbol=Gene) %>% 
	mutate(
		chrom=sapply(strsplit(Variant, split=":"), `[`, 1),
		pos=sapply(strsplit(Variant, split=":"), `[`, 2),
		ref_alt=sapply(strsplit(Variant, split=":"), `[`, 3)) %>% 
	mutate(
		ref=sapply(strsplit(ref_alt, split=">"), `[`, 1),
		alt=sapply(strsplit(ref_alt, split=">"), `[`, 2)
	) %>% select(-Variant, -ref_alt)

# Fill in the '-' with the counts, these should be 0 for each of the alt alleles and the 
# cohort size for the reference count.
get_cohort_counts <- function(column, dt) {
	max(names(table(dt[[column]])))
}

melt_dt <- function(dt, val.name) {
	dt <- melt(data.table(dt),
		measure.vars = c(
		"bsc_kaiser_permanente_EUR",
		"bsc_kaiser_permanente_AFR",
		"bsc_kaiser_permanente_LAT",
		"bsc_kaiser_permanente_EAS",
		"bsc_rarebliss",
		"bsc_bridges",
		"bsc_swedish"),
		variable.name="location",
		value.name=val.name)
	return(dt)
}

for (column in c(
	"KP_EUR_Ca_R", "KP_EUR_Co_R",
	"KP_AFR_Ca_R", "KP_AFR_Co_R",
	"KP_LAT_Ca_R", "KP_LAT_Co_R",
	"KP_EAS_Ca_R", "KP_EAS_Co_R",
	"RB_Ca_R", "RB_Co_R",
	"BR_Ca_R", "BR_Co_R",
	"Sw_Ca_R", "Sw_Co_R")
) {
	dt[[gsub("_R$", "_A", column)]][which(dt[[column]] == "-")] <- 0
	dt[[column]][which(dt[[column]] == "-")] <- get_cohort_counts(column, dt)
}

# Select the case columns and melt across locations after renaming columns
dt_case_ref_counts <- dt %>% 
	select(
		chrom, pos, ref, alt, gene_symbol, consequence_category,
		KP_EUR_Ca_R, KP_AFR_Ca_R, KP_LAT_Ca_R, KP_EAS_Ca_R,
		RB_Ca_R, BR_Ca_R, Sw_Ca_R) %>% 
	rename(
		bsc_kaiser_permanente_EUR=KP_EUR_Ca_R,
		bsc_kaiser_permanente_AFR=KP_AFR_Ca_R,
		bsc_kaiser_permanente_LAT=KP_LAT_Ca_R,
		bsc_kaiser_permanente_EAS=KP_EAS_Ca_R,
		bsc_rarebliss=RB_Ca_R,
		bsc_bridges=BR_Ca_R,
		bsc_swedish=Sw_Ca_R)

dt_control_ref_counts <- dt %>% 
	select(
		chrom, pos, ref, alt, gene_symbol, consequence_category,
		KP_EUR_Co_R, KP_AFR_Co_R, KP_LAT_Co_R, KP_EAS_Co_R,
		RB_Co_R, BR_Co_R, Sw_Co_R) %>% 
	rename(
		bsc_kaiser_permanente_EUR=KP_EUR_Co_R,
		bsc_kaiser_permanente_AFR=KP_AFR_Co_R,
		bsc_kaiser_permanente_LAT=KP_LAT_Co_R,
		bsc_kaiser_permanente_EAS=KP_EAS_Co_R,
		bsc_rarebliss=RB_Co_R,
		bsc_bridges=BR_Co_R,
		bsc_swedish=Sw_Co_R)

dt_case_alt_counts <- dt %>% 
	select(
		chrom, pos, ref, alt, gene_symbol, consequence_category,
		KP_EUR_Ca_A, KP_AFR_Ca_A, KP_LAT_Ca_A, KP_EAS_Ca_A,
		RB_Ca_A, BR_Ca_A, Sw_Ca_A) %>% 
	rename(
		bsc_kaiser_permanente_EUR=KP_EUR_Ca_A,
		bsc_kaiser_permanente_AFR=KP_AFR_Ca_A,
		bsc_kaiser_permanente_LAT=KP_LAT_Ca_A,
		bsc_kaiser_permanente_EAS=KP_EAS_Ca_A,
		bsc_rarebliss=RB_Ca_A,
		bsc_bridges=BR_Ca_A,
		bsc_swedish=Sw_Ca_A)

dt_control_alt_counts <- dt %>% 
	select(
		chrom, pos, ref, alt, gene_symbol, consequence_category,
		KP_EUR_Co_A, KP_AFR_Co_A, KP_LAT_Co_A, KP_EAS_Co_A,
		RB_Co_A, BR_Co_A, Sw_Co_A) %>% 
	rename(
		bsc_kaiser_permanente_EUR=KP_EUR_Co_A,
		bsc_kaiser_permanente_AFR=KP_AFR_Co_A,
		bsc_kaiser_permanente_LAT=KP_LAT_Co_A,
		bsc_kaiser_permanente_EAS=KP_EAS_Co_A,
		bsc_rarebliss=RB_Co_A,
		bsc_bridges=BR_Co_A,
		bsc_swedish=Sw_Co_A)

dt_case_ref_counts <- melt_dt(dt_case_ref_counts, "burden.is_BP.case_no_count")
dt_control_ref_counts <-  melt_dt(dt_control_ref_counts, "burden.is_BP.control_no_count")
dt_case_alt_counts <- melt_dt(dt_case_alt_counts, "burden.is_BP.case_count")
dt_control_alt_counts <- melt_dt(dt_control_alt_counts, "burden.is_BP.control_count")

dt_case <- merge(dt_case_ref_counts, dt_case_alt_counts)
dt_control <- merge(dt_control_ref_counts, dt_control_alt_counts)
dt <- merge(dt_case, dt_control)

fwrite(dt, file="../../BSC_data/MAC5/parsed_and_combined_bsc_MAC5.tsv", sep='\t')
# Move to the cloud
system("gsutil cp ../../BSC_data/MAC5/parsed_and_combined_bsc_MAC5.tsv gs://raw_data_bipolar_dalio_w1_w2/inputs/BSC_MAC5_counts.tsv")

