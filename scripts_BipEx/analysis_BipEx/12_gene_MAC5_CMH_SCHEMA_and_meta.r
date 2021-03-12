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

# This time round, given that we have very diverse ancestries, we shouldn't run a Fisher's Exact test.
# Just run a CMH for each of the genes.
# Note also, that we do not have the non_coding category, so need to do the test controlling for overall burden that is not non-coding.

# # Previous version which didn't include UK/Ireland cases and controls
# dt_schema <- fread("../../SCHEMA_data/schema-counts-grch37-no-ukirl.tsv", header=TRUE) %>% 
# 	filter(csq=='lof') %>% 
# 	dplyr::rename(
# 		consequence_category = csq, is_SCZ.case_count = `TRUE`, is_SCZ.control_count = `FALSE`,
# 		gene_symbol = gene_name,
# 		forest_location = group)

# Latest version that doesn't include UK/Ireland controls
dt_schema <- fread("../../SCHEMA_data/schema-counts-grch37-no-ukirl-controls.tsv", header=TRUE) %>% 
	filter(csq=='lof') %>% 
	dplyr::rename(
		consequence_category = csq, is_SCZ.case_count = `TRUE`, is_SCZ.control_count = `FALSE`,
		gene_symbol = gene_name,
		forest_location = group)

# Add in the case and control 'no count' data.

# # Previous version which didn't include UK/Ireland cases and controls
# dt_strata <- fread("../../SCHEMA_data/schema-no-ukirl-strata-counts.tsv", header=TRUE) %>% rename(forest_location = group, case_count=`1`, control_count=`0`)
# Latest version that doesn't include UK/Ireland controls
dt_strata <- fread("../../SCHEMA_data/schema-strata-counts-no-ukirl-controls.tsv", header=TRUE) %>% dplyr::rename(forest_location = group, case_count=`1`, control_count=`0`)

dt_schema <- merge(dt_schema, dt_strata, by='forest_location') %>% 
	mutate(
		is_SCZ.case_no_count = case_count - is_SCZ.case_count,
		is_SCZ.control_no_count = control_count - is_SCZ.control_count
		) %>% select(-case_count, -control_count)

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
			dt_obs <- dt_obs %>% dplyr::rename(pval = pval_cc, OR = R)
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


# This should be Fisher's exact, plus Stouffer on CMH for SCHEMA.
# As it stands, we are not as well powered as we could be.
# Evaluate the p-values using CMH in SCHEMA and then meta-analyse with our Fisher's exact tests.

# First get the counts for BP before adding in the data and evaluate CMH and Fisher's exact tests.
consequences <- c("lof")

# For CMH tests - MAC5
dt_to_merge <- dt_schema %>% filter(consequence_category=="lof") %>% select(gene_symbol, gene_id)
# Filter to the unique rows
dt_to_merge <- unique(dt_to_merge)

dt_gene_test_schema <- create_CMH_supplementary_table(dt_schema, "is_SCZ") %>% 
	select(gene_symbol, consequence_category, is_SCZ.pval, is_SCZ.OR)
dt_gene_test_schema <- merge(dt_gene_test_schema, dt_to_merge, by="gene_symbol")
# I have the gene-test results for BipEx.
# Read in the results and meta-analyse the p-values where the two are available.

load("../../analysis_plots/gene_counts_qq/plots/BP_gene_fisher_MAC5_gnom_non_psych_qq_gnom_non_psych.is_BP.Rdata")
# Read in the large count data file so that I can merge on ENSG gene_id
dt_to_merge <- fread('gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/14_gene_symbol_gene_id_key.tsv.bgz | zcat')
dt_gene_test_bipex <- dt_plot %>% mutate(pval=10^(-pval)) %>% filter(consequence_category == "ptv") %>% dplyr::rename(gene_symbol=labels)
dt_gene_test_bipex <- merge(dt_gene_test_bipex, dt_to_merge, by='gene_symbol')

dt_gene_test <- merge(dt_gene_test_bipex, dt_gene_test_schema, by='gene_id') %>% 
dplyr::rename(gene_symbol=gene_symbol.x, consequence_category=consequence_category.x) %>% select(-gene_symbol.y, -consequence_category.y)

# Now meta-analyse using Stouffer's method.
# We shouldn't actually look at the others, just really look at the collection of genes that are at the top of the BP1 list.

# Sanity check that before is correct.
dt_gene_test <- dt_gene_test %>% mutate(meta_stouffer_p = pnorm((qnorm(pval) + qnorm(is_SCZ.pval))/sqrt(2)))
dt_gene_test <- dt_gene_test %>% arrange(pval)

# Now add in the case and control counts (summed across strata in the case of SCHEMA).
dt_schema_counts <- dt_schema %>% filter(consequence_category=="lof") %>% 
dplyr::group_by(gene_id) %>% 
dplyr::summarize(
	SCHEMA_case_no_count=sum(is_SCZ.case_no_count),
	SCHEMA_case_count=sum(is_SCZ.case_count),
	SCHEMA_control_no_count=sum(is_SCZ.control_no_count),
	SCHEMA_control_count=sum(is_SCZ.control_count)) %>% select(
	gene_id,
	SCHEMA_case_no_count, SCHEMA_case_count,
	SCHEMA_control_no_count, SCHEMA_control_count)

# Get the BipEx counts to merge in.
GENE_OUT_BP_including_BPSCZ_MAC5_TSV <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_bool.tsv.bgz | gzcat'
dt_bipex_counts <- fread(GENE_OUT_BP_including_BPSCZ_MAC5_TSV) %>% 
dplyr::filter(consequence_category == "ptv") %>% 
dplyr::group_by(gene_symbol) %>% 
dplyr::summarize(
	BipEx_BD_case_no_count =sum(burden_gnom_non_psych.is_BP.case_no_count),
	BipEx_BD_case_count =sum(burden_gnom_non_psych.is_BP.case_count),
	BipEx_control_no_count =sum(burden_gnom_non_psych.is_BP.control_no_count),
	BipEx_control_count =sum(burden_gnom_non_psych.is_BP.control_count),
	BipEx_BD1_case_no_count =sum(burden_gnom_non_psych.is_BP1.case_no_count),
	BipEx_BD1_case_count =sum(burden_gnom_non_psych.is_BP1.case_count),
	BipEx_control_no_count =sum(burden_gnom_non_psych.is_BP1.control_no_count),
	BipEx_control_count =sum(burden_gnom_non_psych.is_BP1.control_count),
	BipEx_BD2_case_no_count =sum(burden_gnom_non_psych.is_BP2.case_no_count),
	BipEx_BD2_case_count =sum(burden_gnom_non_psych.is_BP2.case_count),
	BipEx_control_no_count =sum(burden_gnom_non_psych.is_BP2.control_no_count),
	BipEx_control_count =sum(burden_gnom_non_psych.is_BP2.control_count),
	) %>% select(
	gene_symbol,
	BipEx_BD_case_no_count,
	BipEx_BD_case_count,
	BipEx_control_no_count,
	BipEx_control_count,
	BipEx_BD1_case_no_count,
	BipEx_BD1_case_count,
	BipEx_control_no_count,
	BipEx_control_count,
	BipEx_BD2_case_no_count,
	BipEx_BD2_case_count,
	BipEx_control_no_count,
	BipEx_control_count)

dt_gene_test <- merge(merge(dt_gene_test, dt_bipex_counts, by='gene_symbol'), dt_schema_counts, by='gene_id') %>% 
mutate(case_count = SCHEMA_case_count + BipEx_BD_case_count, control_count = SCHEMA_control_count + BipEx_control_count, 
	   case_no_count = SCHEMA_case_no_count + BipEx_BD_case_no_count, control_no_count =  SCHEMA_control_no_count + BipEx_control_no_count) %>% 
mutate(`combined OR` = (case_count / control_count) / (case_no_count / control_no_count),
	SCHEMA_case = SCHEMA_case_count + SCHEMA_case_no_count,
	SCHEMA_control= SCHEMA_control_count + SCHEMA_control_no_count,
	BipEx_BD_case = BipEx_BD_case_count + BipEx_BD_case_no_count,
	BipEx_BD_control = BipEx_control_count + BipEx_control_no_count) %>% 
mutate(SCHEMA_N = SCHEMA_case + SCHEMA_control,
	BipEx_N = BipEx_BD_case + BipEx_BD_control, 
	SCHEMA_Neff = (4 * (SCHEMA_N) * (SCHEMA_case / SCHEMA_N) * (SCHEMA_control / SCHEMA_N)), 
	BipEx_Neff = (4 * (BipEx_N) * (BipEx_BD_case / BipEx_N) * (BipEx_BD_control / BipEx_N)),
	SCHEMA_weight = sqrt(SCHEMA_Neff), BipEx_weight = sqrt(BipEx_Neff), 
	BipEx_Z = sign(OR) * qnorm(pval/2), SCHEMA_Z = sign(is_SCZ.OR) * qnorm(is_SCZ.pval/2),
	weighted_Z = (SCHEMA_weight * SCHEMA_Z + BipEx_weight * BipEx_Z) / sqrt(SCHEMA_weight^2 + BipEx_weight^2),
	meta_p = 2 * pnorm(weighted_Z) 
	) %>% arrange(pval)
dt_gene_test_to_write <- dt_gene_test %>%
select(gene_symbol, gene_id,
	BipEx_BD_case_count, BipEx_BD1_case_count, BipEx_BD2_case_count, BipEx_control_count,
	BipEx_BD_case_no_count, BipEx_BD1_case_no_count, BipEx_BD2_case_no_count, BipEx_control_no_count,
	pval, OR,
	SCHEMA_case_count, SCHEMA_control_count,
	SCHEMA_case_no_count, SCHEMA_control_no_count,
	BipEx_BD_case, BipEx_BD_control,
	SCHEMA_case, SCHEMA_control,
	is_SCZ.pval, is_SCZ.OR, `combined OR`, meta_p) %>%
dplyr::rename(
	Gene=gene_symbol,
	`Meta p-value`=meta_p,
	`BipEx BD p-value`=pval, 
	`BipEx BD OR`=OR,
	`SCHEMA p-value`=is_SCZ.pval, 
	`SCHEMA OR`=is_SCZ.OR,
	) %>% arrange(`BipEx BD p-value`)

fwrite(head(dt_gene_test_to_write, 10), "../../paper_tables/table2.tsv", sep='\t')
fread("../../paper_tables/table2.tsv")


