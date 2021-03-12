library(dplyr)
library(data.table)

load("forest_plot_Rdata_files/BP_including_BPSCZ.Rdata")

forest_plots_dt <- rbind(
	forest_plots$is_BP$logit %>% mutate(phenotype="BP"),
	forest_plots$is_BP1$logit %>% mutate(phenotype="BP1"),
	forest_plots$is_BP2$logit %>% mutate(phenotype="BP2")
) %>% filter(grepl('gnom_non_psych', label)) %>% 
	filter(
		grepl('pli_09\\.', label) & 
		(
			grepl('synonymous', label) |
			grepl('MPC_2_missense', label) |
			grepl('other_missense', label) |
			grepl('PTV', label)
		)
	)

fwrite(forest_plots_dt, file='table_for_TJ.tsv', sep='\t')
