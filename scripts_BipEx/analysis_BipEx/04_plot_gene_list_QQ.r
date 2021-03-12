library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrastr)

dt <- fread("gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/firth_coding_burden_BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.tsv.bgz | gzcat")

phenotypes <- unique(dt$phenotype)
log10Pe <- expression(paste("Expected -log"[10], "(", italic(p), ")"))
log10Po <- expression(paste("Observed -log"[10], "(", italic(p), ")"))

ci <- 0.95

dt <- dt %>% mutate(p_value = ifelse(beta > 0, p_value/2, 1-(p_value)/2))
dt <- dt %>% filter(!(gene_set %in% unique(grep("constrain", gene_set, value=TRUE))))
dt$consequence_category <- mapvalues(dt$consequence_category,
	from = c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv"),
	to = c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV"))
dt <- dt %>% filter(consequence_category != "Non-coding")
dt$consequence_category <- ordered(dt$consequence_category, levels = c("Synonymous", "Other missense", "Damaging missense", "PTV"))

pdf(file="QQ_genelist_syn.pdf", width=4, height=4)
for (phe in phenotypes) {
	dt_current <- dt %>% filter(phenotype == phe) %>% mutate(log_p=-log10(p_value))
	dt_current <- data.table(dt_current)
	setkeyv(dt_current, c("consequence_category", "p_value"))
	n_consequence_categories <- length(unique(dt_current$consequence_category))
	n_genesets <- nrow(dt_current) / n_consequence_categories

	dt_current[, expected:=rep(-log10(seq(1, n_genesets)/(n_genesets + 1)), n_consequence_categories)]
	dt_current[, clower:=rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]
	dt_current[, cupper:=rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]

	p <- ggplot(dt_current, aes(x=expected,y=log_p, color=consequence_category)) +
		geom_ribbon(aes(ymin=clower, ymax=cupper), fill = "grey80", color="grey80") +
		geom_point_rast() + theme_classic() + xlab(log10Pe) + ylab(log10Po) +
		scale_color_manual(values=c("#FCBBA1", "#FC9272", "#FB6A4A", "#CB181D")) + 
		theme(legend.title=element_blank(), legend.position="none")  + geom_abline(intercept=0, slope=1, size=1)
	print(p)
}
dev.off()

dt <- fread("gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/firth_coding_burden_BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.tsv.bgz | gzcat")

phenotypes <- unique(dt$phenotype)
log10Pe <- expression(paste("Expected -log"[10], "(", italic(p), ")"))
log10Po <- expression(paste("Observed -log"[10], "(", italic(p), ")"))

ci <- 0.95

dt <- dt %>% filter(!(gene_set %in% unique(
		c(grep("constrain", gene_set, value=TRUE),
		grep("KEGG", gene_set, value=TRUE),
		grep("PANTHER", gene_set, value=TRUE),
		grep("REACTOME", gene_set, value=TRUE),
		grep("GOMF", gene_set, value=TRUE),
		grep("GOCC", gene_set, value=TRUE),
		grep("GOBP", gene_set, value=TRUE),
		grep("HALLMARK", gene_set, value=TRUE))
	)))
dt <- dt %>% filter(!(gene_set %in% unique(grep("constrain", gene_set, value=TRUE))))
dt$consequence_category <- mapvalues(dt$consequence_category,
	from = c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv"),
	to = c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV"))
dt <- dt %>% filter(consequence_category != "Non-coding")
dt$consequence_category <- ordered(dt$consequence_category, levels = c("Synonymous", "Other missense", "Damaging missense", "PTV"))

pdf(file="QQ_genelist_syn_targeted.pdf", width=4, height=4)
for (phe in phenotypes) {
	dt_current <- dt %>% filter(phenotype == phe) %>% mutate(log_p=-log10(p_value))
	dt_current <- data.table(dt_current)
	setkeyv(dt_current, c("consequence_category", "p_value"))
	n_consequence_categories <- length(unique(dt_current$consequence_category))
	n_genesets <- nrow(dt_current) / n_consequence_categories

	dt_current[, expected:=rep(-log10(seq(1, n_genesets)/(n_genesets + 1)), n_consequence_categories)]
	dt_current[, clower:=rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]
	dt_current[, cupper:=rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]

	p <- ggplot(dt_current, aes(x=expected,y=log_p, color=consequence_category)) +
		geom_ribbon(aes(ymin=clower, ymax=cupper), fill = "grey80", color="grey80") +
		geom_point_rast() + theme_classic() + xlab(log10Pe) + ylab(log10Po) +
		scale_color_manual(values=c("#FCBBA1", "#FC9272", "#FB6A4A", "#CB181D")) + 
		theme(legend.title=element_blank(), legend.position="none")  + geom_abline(intercept=0, slope=1, size=1)
	print(p)
}
dev.off()

dt <- fread("gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/firth_coding_burden_BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.tsv.bgz | gzcat")

phenotypes <- unique(dt$phenotype)
log10Pe <- expression(paste("Expected -log"[10], "(", italic(p), ")"))
log10Po <- expression(paste("Observed -log"[10], "(", italic(p), ")"))

ci <- 0.95

dt <- dt %>% mutate(p_value = ifelse(beta > 0, p_value/2, 1-(p_value)/2))
dt <- dt %>% filter(gene_set %in% unique(
	c(grep("KEGG", gene_set, value=TRUE),
	  grep("PANTHER", gene_set, value=TRUE),
	  grep("REACTOME", gene_set, value=TRUE),
	  grep("GOMF", gene_set, value=TRUE),
	  grep("GOCC", gene_set, value=TRUE),
	  grep("GOBP", gene_set, value=TRUE),
	  grep("HALLMARK", gene_set, value=TRUE))
	))
dt$consequence_category <- mapvalues(dt$consequence_category,
	from = c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv"),
	to = c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV"))
dt <- dt %>% filter(consequence_category != "Non-coding")
dt$consequence_category <- ordered(dt$consequence_category, levels = c("Synonymous", "Other missense", "Damaging missense", "PTV"))

pdf(file="QQ_genelist_syn_non_targeted.pdf", width=4, height=4)
for (phe in phenotypes) {
	dt_current <- dt %>% filter(phenotype == phe) %>% mutate(log_p=-log10(p_value))
	dt_current <- data.table(dt_current)
	setkeyv(dt_current, c("consequence_category", "p_value"))
	n_consequence_categories <- length(unique(dt_current$consequence_category))
	n_genesets <- nrow(dt_current) / n_consequence_categories

	dt_current[, expected:=rep(-log10(seq(1, n_genesets)/(n_genesets + 1)), n_consequence_categories)]
	dt_current[, clower:=rep(-log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]
	dt_current[, cupper:=rep(-log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n_genesets, shape2 = n_genesets:1)), n_consequence_categories)]

	p <- ggplot(dt_current, aes(x=expected,y=log_p, color=consequence_category)) +
		geom_ribbon(aes(ymin=clower, ymax=cupper), fill = "grey80", color="grey80") +
		geom_point_rast() + theme_classic() + xlab(log10Pe) + ylab(log10Po) +
		scale_color_manual(values=c("#FCBBA1", "#FC9272", "#FB6A4A", "#CB181D")) + 
		theme(legend.title=element_blank(), legend.position="none") + geom_abline(intercept=0, slope=1, size=1)
	print(p)
}
dev.off()

# Finally, extract the p-values for the SCZ gene-lists that we just passed.

dt <- fread("gsutil cat gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/firth_coding_burden_BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists_SCZ_overlap.tsv.bgz | gzcat")

phenotypes <- unique(dt$phenotype)
dt$consequence_category <- mapvalues(dt$consequence_category,
	from = c("non_coding", "synonymous", "other_missense", "damaging_missense", "ptv"),
	to = c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV"))
dt <- dt %>% filter(consequence_category != "Non-coding")
dt$consequence_category <- ordered(dt$consequence_category, levels = c("Synonymous", "Other missense", "Damaging missense", "PTV"))

