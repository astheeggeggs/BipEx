library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggrepel)
library(rjson)

source('r_functions/burden_tests.r')

# First, we need to parse the data and put it in the following form:

# gene_symbol, consequence_category, 
# burden.is_BP.case_count, burden.is_BP.case_no_count, burden.is_BP.control_count, burden.is_BP.control_no_count,
# burden_gnom_non_psych.is_BP.case_count, burden_gnom_non_psych.is_BP.case_no_count, burden_gnom_non_psych.is_BP.control_count, burden_gnom_non_psych.is_BP.control_no_count, 

# For the Fisher's exact tests - two versions - one with MAC5, and one singleton version.
# and the same but also including location (for us this will be the locations provided by the BSC).

# Note also, that we do not have the non_coding category, so need to do the test controlling for overall burden that is not non-coding.

GENE_OUT_BP_including_BPSCZ_MAC5_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool.tsv.bgz | gzcat'
GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool_by_location.tsv.bgz | gzcat'

dt_bsc <- fread("gsutil cat gs://raw_data_bipolar_dalio_w1_w2/analysis/03_BSC_MAC5_GRCh38_and_gnomAD.tsv")

dt_bsc_location <- dt_bsc %>% 
select(-chrom, -pos, -ref, -alt, -locus, -alleles, -new_locus.result, -new_locus.is_negative_strand,
	-new_alleles) %>% rename(forest_location=location)
# Ensure no gene_symbol changes
dt_bipex <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_TSV)
all(dt_bsc_location$gene_symbol %in% dt_bipex$gene_symbol)
# Yes, good.

# Group and sum.
dt_bsc_location <- dt_bsc_location %>% mutate(
	case_count=burden.is_BP.case_count + burden.is_BP.case_no_count,
	control_count=burden.is_BP.control_count + burden.is_BP.control_no_count,
	burden_gnom_non_psych.is_BP.case_count= burden.is_BP.case_count * !inGnomAD_nonpsych,
	burden_gnom_non_psych.is_BP.control_count= burden.is_BP.control_count * !inGnomAD_nonpsych
)

dt_bsc_location$consequence_category[which(dt_bsc_location$consequence_category == "damaging")] <- "damaging_missense"

# Filter down to get the collection of variants that are MAC <= 5 in the BSC dataset.
dt_bsc_gene_burden_MAC5 <- dt_bsc_location %>% 
	group_by(gene_symbol, consequence_category, forest_location) %>% 
	summarise(
		burden.is_BP.case_count=sum(burden.is_BP.case_count),
		burden.is_BP.control_count=sum(burden.is_BP.control_count),
		# Now those with the additional restriction that they are not in gnomAD.
		burden_gnom_non_psych.is_BP.case_count=sum(burden_gnom_non_psych.is_BP.case_count),
		burden_gnom_non_psych.is_BP.control_count=sum(burden_gnom_non_psych.is_BP.control_count),
		# Counts of samples in each gene contingency table.
		case_count=max(case_count),
		control_count=max(control_count)
	) %>% mutate(
		burden.is_BP.case_no_count = case_count - burden.is_BP.case_count,
		burden.is_BP.control_no_count = control_count - burden.is_BP.control_count,
		burden_gnom_non_psych.is_BP.case_no_count = case_count - burden_gnom_non_psych.is_BP.case_count,
		burden_gnom_non_psych.is_BP.control_no_count = control_count - burden_gnom_non_psych.is_BP.control_count
	) %>% select(-case_count, -control_count)

# Now write the results.
fwrite(dt_bsc_gene_burden_MAC5 , "../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_MAC5_gene_bool_by_location.tsv")

# Sum across locations and write (this is for the Fisher's exact test)
fwrite(dt_bsc_gene_burden_MAC5 %>% 
	select(-forest_location) %>% 
	group_by(gene_symbol, consequence_category) %>%
	summarise(
		burden.is_BP.case_count=sum(burden.is_BP.case_count),
		burden.is_BP.control_count=sum(burden.is_BP.control_count),
		burden.is_BP.case_no_count=sum(burden.is_BP.case_no_count),
		burden.is_BP.control_no_count=sum(burden.is_BP.control_no_count),
		burden_gnom_non_psych.is_BP.case_count=sum(burden_gnom_non_psych.is_BP.case_count),
		burden_gnom_non_psych.is_BP.control_count=sum(burden_gnom_non_psych.is_BP.control_count),
		burden_gnom_non_psych.is_BP.case_no_count=sum(burden_gnom_non_psych.is_BP.case_no_count),
		burden_gnom_non_psych.is_BP.control_no_count=sum(burden_gnom_non_psych.is_BP.control_no_count)
		) , "../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_MAC5_gene_bool.tsv")

create_fisher_supplementary_table <- function(dt, phenotype_tests, log=FALSE)
{
	init <- TRUE
	for (phenotype in phenotype_tests)
	{	
		cat(paste0(phenotype, "..."))
		cat("obtaining observed p-values...\n")
		dt_obs <- get_fisher_genes(phenotype, dt) %>% 
		select(gene_symbol, consequence_category, pval, OR)

		if (log == TRUE) {
			dt_obs <- dt_obs %>% mutate(pval = log10(pval), log_OR = log10(OR)) %>% select(-OR)
			names_to_append <- which(names(dt_obs) %in% c("pval", "log_OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("log_pval", "log_OR"))
		} else {
			names_to_append <- which(names(dt_obs) %in% c("pval", "OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("pval", "OR"))
		}
		cat("obtained observed p-values...\n")
		
		if (init) {
			dt_full <- merge(dt, dt_obs, by=c("gene_symbol", "consequence_category"))
			init <- FALSE
		} else {
			dt_full <- merge(dt_full, dt_obs, by=c("gene_symbol", "consequence_category"))
		}
	}
	return(data.table(dt_full))
}

create_CMH_supplementary_table <- function(dt, phenotype_tests, log=FALSE)
{
	init <- TRUE
	for (phenotype in phenotype_tests)
	{	
		cat(paste0(phenotype, "..."))
		cat("obtaining observed p-values...\n")
		dt_obs <- get_CMH_genes(phenotype, dt) %>% 
		select(gene_symbol, consequence_category, pval_cc, R)

		if (log == TRUE) {
			dt_obs <- dt_obs %>% 
				mutate(pval = -log10(pval_cc), log_OR = log10(R)) %>% 
				select(-R, -pval_cc)
			names_to_append <- which(names(dt_obs) %in% c("pval", "log_OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("log_pval", "log_OR"))
		} else {
			dt_obs <- dt_obs %>% rename(pval = pval_cc, OR = R)
			names_to_append <- which(names(dt_obs) %in% c("pval", "OR"))
			names(dt_obs)[names_to_append] <- paste0(phenotype, ".", c("pval", "OR"))
		}
		cat("obtained observed p-values...\n")
		
		if (init) {
			dt_full <- dt_obs
			init <- FALSE
		} else {
			dt_full <- merge(dt_full, dt_obs, by=c("gene_symbol", "consequence_category"))
		}
	}
	return(data.table(dt_full))
}


# Just do a before and after test.

# Note - the collection of genes that we asked for were the top genes for SINGLETONS without the gnomAD restriction.
# we can tell this from the QQ plot on the left here: https://astheeggeggs.github.io/BipEx/gene_analyses_singleton.html,
# and comparing to the list of genes in the variant list.

# First get the counts for BP before adding in the data and evaluate CMH and Fisher's exact tests.
consequences <- c("synonymous", "other_missense", "damaging_missense", "ptv")

# For Fisher's exact tests
dt_bsc_MAC5_fisher <- fread("../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_MAC5_gene_bool.tsv")

# For CMH tests
dt_bsc_MAC5_CMH <- fread("../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_MAC5_gene_bool_by_location.tsv")

# Read in and clean up
dt_bipex_MAC5 <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_TSV) %>% 
	filter(consequence_category %in% consequences)
names(dt_bipex_MAC5) <- gsub("burden.", "", names(dt_bipex_MAC5)) 
names(dt_bipex_MAC5) <- gsub("burden_", "", names(dt_bipex_MAC5))

dt_bipex_MAC5 <- dt_bipex_MAC5 %>% 
	select(gene_symbol, consequence_category,
		is_BP.case_count, is_BP.case_no_count, is_BP.control_count, is_BP.control_no_count,
		gnom_non_psych.is_BP.case_count, gnom_non_psych.is_BP.case_no_count,
		gnom_non_psych.is_BP.control_count, gnom_non_psych.is_BP.control_no_count
		)

dt_bsc_bipex_MAC5 <- merge(dt_bsc_MAC5_fisher, dt_bipex_MAC5) %>%
mutate(
	is_BP.case_count = burden.is_BP.case_count + is_BP.case_count,
	is_BP.control_count = burden.is_BP.control_count + is_BP.control_count,
	is_BP.case_no_count = burden.is_BP.case_no_count + is_BP.case_no_count,
	is_BP.control_no_count = burden.is_BP.control_no_count + is_BP.control_no_count,

	gnom_non_psych.is_BP.case_count = burden_gnom_non_psych.is_BP.case_count + gnom_non_psych.is_BP.case_count,
	gnom_non_psych.is_BP.control_count = burden_gnom_non_psych.is_BP.control_count + gnom_non_psych.is_BP.control_count,
	gnom_non_psych.is_BP.case_no_count = burden_gnom_non_psych.is_BP.case_no_count + gnom_non_psych.is_BP.case_no_count,
	gnom_non_psych.is_BP.control_no_count = burden_gnom_non_psych.is_BP.control_no_count + gnom_non_psych.is_BP.control_no_count
	) %>% select(
	-burden.is_BP.case_count, -burden.is_BP.control_count, -burden.is_BP.case_no_count, -burden.is_BP.control_no_count,
	-burden_gnom_non_psych.is_BP.case_count, -burden_gnom_non_psych.is_BP.control_count,
	-burden_gnom_non_psych.is_BP.case_no_count, -burden_gnom_non_psych.is_BP.control_no_count
	)

# Then merge in the counts from the BSC and rerun.

# MAC5
dt_before_MAC5 <- merge(
	create_fisher_supplementary_table(dt_bipex_MAC5, "is_BP") %>% select(
		gene_symbol, consequence_category,
		is_BP.case_count, is_BP.control_count, is_BP.pval, is_BP.OR),
	create_fisher_supplementary_table(dt_bipex_MAC5, "gnom_non_psych.is_BP") %>% select(
		gene_symbol, consequence_category,
		gnom_non_psych.is_BP.case_count, gnom_non_psych.is_BP.control_count, gnom_non_psych.is_BP.pval, gnom_non_psych.is_BP.OR), by=c("gene_symbol", "consequence_category")
	)

# MAC5
dt_after_MAC5 <- merge(
	create_fisher_supplementary_table(dt_bsc_bipex_MAC5, "is_BP") %>% select(
		gene_symbol, consequence_category,
		is_BP.case_count, is_BP.control_count, is_BP.pval, is_BP.OR),
	create_fisher_supplementary_table(dt_bsc_bipex_MAC5, "gnom_non_psych.is_BP") %>% select(
		gene_symbol, consequence_category,
		gnom_non_psych.is_BP.case_count, gnom_non_psych.is_BP.control_count, gnom_non_psych.is_BP.pval, gnom_non_psych.is_BP.OR), by=c("gene_symbol", "consequence_category")
	)

# Now, merge these all together.
dt_before_and_after <- merge(dt_before, dt_after, by=c("gene_symbol", "consequence_category"))
names(dt_before_and_after) <- gsub("\\.x", ".before", names(dt_before_and_after))
names(dt_before_and_after) <- gsub("\\.y", ".after", names(dt_before_and_after))
dt_before_and_after_MAC5 <- merge(dt_before_MAC5, dt_after_MAC5, by=c("gene_symbol", "consequence_category"))
names(dt_before_and_after_MAC5) <- gsub("\\.x", ".before", names(dt_before_and_after_MAC5))
names(dt_before_and_after_MAC5) <- gsub("\\.y", ".after", names(dt_before_and_after_MAC5))

fwrite(dt_before_and_after_MAC5, file="gene_tables/before_and_after_Fisher_bsc.tsv", sep='\t')

# Do the same for the CMH using all the BSC locations as strata.
consequences <- c("synonymous", "other_missense", "damaging_missense", "ptv")

dt_bsc_MAC5_CMH <- fread("../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_MAC5_gene_bool_by_location.tsv")
dt_bsc_CMH <- fread("../../BSC_data/burden_output_for_fisher_and_CMH/03_BSC_BP_singleton_gene_bool_by_location.tsv")

# Read in and clean up
dt_bipex_MAC5 <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_BY_LOCATION_TSV) %>% 
	filter(consequence_category %in% consequences)
names(dt_bipex_MAC5) <- gsub("burden.", "", names(dt_bipex_MAC5)) 
names(dt_bipex_MAC5) <- gsub("burden_", "", names(dt_bipex_MAC5))

dt_bipex_MAC5 <- dt_bipex_MAC5 %>% 
	select(gene_symbol, consequence_category, forest_location,
		is_BP.case_count, is_BP.case_no_count, is_BP.control_count, is_BP.control_no_count,
		gnom_non_psych.is_BP.case_count, gnom_non_psych.is_BP.case_no_count,
		gnom_non_psych.is_BP.control_count, gnom_non_psych.is_BP.control_no_count
		)

dt_bsc_bipex_MAC5 <- merge(dt_bsc_MAC5_CMH %>% 
	rename(
		is_BP.case_count = burden.is_BP.case_count,
		is_BP.control_count = burden.is_BP.control_count,
		is_BP.case_no_count = burden.is_BP.case_no_count,
		is_BP.control_no_count = burden.is_BP.control_no_count,
		gnom_non_psych.is_BP.case_count = burden_gnom_non_psych.is_BP.case_count,
		gnom_non_psych.is_BP.control_count = burden_gnom_non_psych.is_BP.control_count,
		gnom_non_psych.is_BP.case_no_count = burden_gnom_non_psych.is_BP.case_no_count,
		gnom_non_psych.is_BP.control_no_count = burden_gnom_non_psych.is_BP.control_no_count
		), dt_bipex_MAC5, all=TRUE)

dt_before_MAC5 <- merge(
	create_CMH_supplementary_table(dt_bipex_MAC5, "is_BP") %>% select(
		gene_symbol, consequence_category, is_BP.pval, is_BP.OR),
	create_CMH_supplementary_table(dt_bipex_MAC5, "gnom_non_psych.is_BP") %>% select(
		gene_symbol, consequence_category, gnom_non_psych.is_BP.pval, gnom_non_psych.is_BP.OR), by=c("gene_symbol", "consequence_category")
	)

dt_after_MAC5 <- merge(
	create_CMH_supplementary_table(dt_bsc_bipex_MAC5, "is_BP") %>% select(
		gene_symbol, consequence_category, is_BP.pval, is_BP.OR),
	create_CMH_supplementary_table(dt_bsc_bipex_MAC5, "gnom_non_psych.is_BP") %>% select(
		gene_symbol, consequence_category, gnom_non_psych.is_BP.pval, gnom_non_psych.is_BP.OR), by=c("gene_symbol", "consequence_category"),
	all.x=TRUE)

# Now, merge these all together.
dt_before_and_after_MAC5 <- merge(dt_before_MAC5, dt_after_MAC5, by=c("gene_symbol", "consequence_category"))
names(dt_before_and_after_MAC5) <- gsub("\\.x", ".before", names(dt_before_and_after_MAC5))
names(dt_before_and_after_MAC5) <- gsub("\\.y", ".after", names(dt_before_and_after_MAC5))

# Get the BipEx counts to merge in.
GENE_OUT_BP_including_BPSCZ_MAC5_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool.tsv.bgz | gzcat'
dt_bipex_counts <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_TSV) %>% 
dplyr::filter(consequence_category == "ptv") %>% 
dplyr::group_by(gene_symbol) %>% 
dplyr::summarize(
	`BipEx BD case no count`=sum(burden_gnom_non_psych.is_BP.case_no_count),
	`BipEx BD case count`=sum(burden_gnom_non_psych.is_BP.case_count),
	`BipEx control no count`=sum(burden_gnom_non_psych.is_BP.control_no_count),
	`BipEx control count`=sum(burden_gnom_non_psych.is_BP.control_count)
	) %>% select(
	gene_symbol,
	`BipEx BD case count`,
	`BipEx control count`)

# Get the BSC counts to merge in.
dt_bsc_counts <- dt_bsc_MAC5_CMH %>% filter(consequence_category=="ptv") %>%
dplyr::group_by(gene_symbol) %>% 
dplyr::summarize(
	`BSC BD case no count`=sum(burden_gnom_non_psych.is_BP.case_no_count),
	`BSC BD case count`=sum(burden_gnom_non_psych.is_BP.case_count),
	`BSC control no count`=sum(burden_gnom_non_psych.is_BP.control_no_count),
	`BSC control count`=sum(burden_gnom_non_psych.is_BP.control_count)
	) %>% select(
	gene_symbol,
	`BSC BD case count`,
	`BSC control count`)

# Just write the not in gnomAD data to a tsv.
fwrite(dt_before_and_after_MAC5, file="gene_tables/before_and_after_bsc.tsv", sep='\t')
table1 <- dt_before_and_after_MAC5 %>% 
	filter(consequence_category=="ptv") %>% 
	select(
		gene_symbol,
		gnom_non_psych.is_BP.pval.before,
		gnom_non_psych.is_BP.OR.before,
		gnom_non_psych.is_BP.pval.after,
		gnom_non_psych.is_BP.OR.after
		) %>%
	rename(
		`BipEx BP p-value`=gnom_non_psych.is_BP.pval.before,
		`BipEx BP OR`=gnom_non_psych.is_BP.OR.before,
		`BipEx & BSC BP p-value`=gnom_non_psych.is_BP.pval.after,
		`BipEx & BSC BP OR`=gnom_non_psych.is_BP.OR.after) %>%
	arrange(`BipEx BP p-value`)
table1 <- merge(merge(table1, dt_bipex_counts, by="gene_symbol", all.x=TRUE), dt_bsc_counts, by="gene_symbol", all.x=TRUE) %>% rename(Gene=gene_symbol) %>% arrange(`BipEx BP p-value`)
table1 <- head(table1, 10) %>% select(Gene,
	`BipEx BD case count`,
	`BipEx control count`,
	`BipEx BP p-value`,
	`BipEx BP OR`,
	`BSC BD case count`,
	`BSC control count`,
	`BipEx & BSC BP p-value`,
	`BipEx & BSC BP OR`)
fwrite(table1 , file="../../paper_tables/table1.tsv", sep='\t')
