library(data.table)
library(dplyr)

# Obtain the pLI scores for all of the genes.
dt_pLI <- fread(cmd = "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/14_gene_symbol_gene_id_key.tsv.bgz | gzcat")
dt <- fread("../../analysis_plots/gene_counts_qq/plots/BP_gene_fisher_MAC5_gnom_non_psych_qq.tsv")
dt <- merge(dt, dt_pLI)
head(dt %>% filter(consequence_category=="ptv") %>% filter(pLI > 0.9) %>% arrange(gnom_non_psych.is_BP.pval) %>% select(gene_symbol, gnom_non_psych.is_BP.pval), 20) %>% 
fwrite("../../analysis_plots/gene_counts_qq/plots/top_20_BP_gene_fisher_MAC5_gnom_non_psych_qq.tsv", sep='\t')

system("gsutil cp ../../analysis_plots/gene_counts_qq/plots/top_20_BP_gene_fisher_MAC5_gnom_non_psych_qq.tsv gs://dalio_bipolar_w1_w2_hail_02/data/genes/top_20_BP_gene_fisher_MAC5_gnom_non_psych_qq.tsv")
