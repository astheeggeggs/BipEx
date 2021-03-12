library(data.table)

dt <- fread('gstuil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz | gzcat')
dt_sample <- fread('gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz | gzcat')

dt_SCN3A_ptv <- dt[gene_symbol=='SCN3A' & consequence_category=='ptv']
dt_SCN3A_ptv[, gene_symbol:=NULL]
dt_SCN3A_ptv[, consequence_category:=NULL]
dt_out <- transpose(dt_SCN3A_ptv, keep.names='s')
dt_out <- dt_out %>% filter(V1==1)

dt_out <- data.table(dt_out)
setkey(dt_out, "s")

dt_out <- merge(dt_out, dt_sample) %>% select(s, PROJECT_OR_COHORT, LOCATION, INSTITUTION, PI, PHENOTYPE_COARSE, PHENOTYPE_FINE, PSYCHOSIS) %>% fwrite(file='SCN3A_for_Dennis.tsv', sep='\t')
