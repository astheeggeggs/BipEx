library(data.table)
library(dplyr)
library(ggplot2)
library(ggrastr)

# Want the observed p-values for the genes for synonymous burden, then can hijack this code to create
# plots.
source('r_functions/burden_tests.r')

# Fisher
BP_fisher_dt_files <- c(
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous.tsv",
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous.tsv")

# CMH
BP_CMH_dt_files <- c(
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_50_synonymous_by_location.tsv",
	"null_test_plots/03_BP_including_BPSCZ_MAC5_gene_bool_count_20_synonymous_by_location.tsv")

phenotype <-  "gnom_non_psych.is_BP"
log10Pe <- expression(paste("Expected -log"[10], "(", italic(p), ")"))
log10Po <- expression(paste("Observed -log"[10], "(", italic(p), ")"))
ci <- 0.95

pdf(file="QQ_gene_syn_check.pdf", width=4, height=4)

for(BP_dt_file in BP_CMH_dt_files) {
	dt <- fread(BP_dt_file)
	names(dt) <- gsub("burden.", "", names(dt)) 
	dt_obs <- get_CMH_genes(phenotype, dt) %>% rename(labels=gene_symbol) %>% 
		select(labels, consequence_category, pval_cc, R) %>% 
		mutate(
			log_p = -log10(pval_cc),
			log_OR = log10(R)
			)
	n_genes <- nrow(dt_obs)
	dt_obs[, expected:=-log10(seq(1, n_genes)/(n_genes + 1))]
	dt_obs[, clower:=-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_genes, shape2 = n_genes:1))]
	dt_obs[, cupper:=-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_genes, shape2 = n_genes:1))]
	p <- ggplot(dt_obs, aes(x=expected, y=log_p, colour="indianred3")) +
		geom_ribbon(aes(ymin=clower, ymax=cupper), fill = "grey80", color="grey80") +
		geom_point_rast() + theme_classic() + xlab(log10Pe) + ylab(log10Po) +
		theme(legend.title=element_blank(), legend.position="none") + geom_abline(intercept=0, slope=1, size=1)
	print(p)
}

for(BP_dt_file in BP_fisher_dt_files) {
	dt <- fread(BP_dt_file)
	names(dt) <- gsub("burden.", "", names(dt)) 
	dt_obs <- get_fisher_genes(phenotype, dt) %>% rename(labels=gene_symbol) %>% 
		select(labels, consequence_category, pval, OR) %>% 
		mutate(
			log_p = -log10(pval),
			log_OR = log10(OR)
			)
	n_genes <- nrow(dt_obs)
	dt_obs[, expected:=-log10(seq(1, n_genes)/(n_genes + 1))]
	dt_obs[, clower:=-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_genes, shape2 = n_genes:1))]
	dt_obs[, cupper:=-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_genes, shape2 = n_genes:1))]
	p <- ggplot(dt_obs, aes(x=expected, y=log_p, colour="indianred3")) +
		geom_ribbon(aes(ymin=clower, ymax=cupper), fill = "grey80", color="grey80") +
		geom_point_rast() + theme_classic() + xlab(log10Pe) + ylab(log10Po) +
		theme(legend.title=element_blank(), legend.position="none") + geom_abline(intercept=0, slope=1, size=1)
	print(p)
}

dev.off()
 