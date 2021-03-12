library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

source("../QC_BipEx/r_functions_and_parameters/pretty_plotting.r")

# Generating figures similar to those in Giulio's paper.
dt <- fread('gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/variants/17_final_qc.variants.tsv.bgz | gzcat', sep='\t')
dt <- dt %>% 
	mutate(
		ref_AC = as.numeric(gsub("\\[([0-9]+),.*", "\\1", qc.AC)),
		alt_AC = as.numeric(gsub("\\[.*,([0-9]+).*\\]", "\\1", qc.AC))) %>%
	mutate(MAC = ifelse(ref_AC > alt_AC, alt_AC, ref_AC)) %>% 
	filter(MAC < 7 & MAC > 0)

# Jobs - rename and change the colour scheme.
# Figure 1a)
p <- ggplot(data=dt, aes(x=MAC, fill=inGnomAD_nonpsych)) + geom_bar() + theme_classic() + scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Count", fill = "In GnomAD non-psych") + 
	theme(legend.position = c(0.8, 0.9)) + scale_fill_discrete(labels = c("No", "Yes"))

# Figure 1b)
p <- ggplot(data=dt, aes(x=MAC, fill=consequence_category)) + geom_bar() + theme_classic() + scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Count", fill = "Consequence category") + 
	theme(legend.position = c(0.8, 0.8)) +
	scale_fill_brewer(labels = c(
		"Damaging missense",
		"Non-coding",
		"Other missense",
		"PTV",
		"Synonymous",
		"Unknown"), palette="Set3")
print(p)

# Stacked + percent
p <- ggplot(data= dt, aes(x=MAC, fill=consequence_category)) + geom_bar(position=position_fill()) +
	theme_classic() + scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Proportion", fill = "Consequence category") + 
	theme(legend.position = "none") +
	scale_fill_brewer(labels = c(
		"Damaging missense",
		"Non-coding",
		"Other missense",
		"PTV",
		"Synonymous",
		"Unknown"), palette="Set3")
print(p)

quartz()

p <- ggplot(data= (dt %>% filter(!inGnomAD_nonpsych)), aes(x=MAC, fill=consequence_category)) + geom_bar() + theme_classic() + scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Count", fill = "Consequence category") + 
	theme(legend.position = c(0.8, 0.8)) +
	scale_fill_brewer(labels = c(
		"Damaging missense",
		"Non-coding",
		"Other missense",
		"PTV",
		"Synonymous",
		"Unknown"), palette="Set3")
print(p)

# Stacked + percent
p <- ggplot(data= (dt %>% filter(!inGnomAD_nonpsych)), aes(x=MAC, fill=consequence_category)) + geom_bar(position=position_fill()) +
	theme_classic() + scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Proportion", fill = "Consequence category") + 
	theme(legend.position = "none") +
	scale_fill_brewer(labels = c(
		"Damaging missense",
		"Non-coding",
		"Other missense",
		"PTV",
		"Synonymous",
		"Unknown"), palette="Set3")
print(p)

# Stacked + percent
p <- ggplot(data= (dt %>% filter(inGnomAD_nonpsych)), aes(x=MAC, fill=consequence_category)) + geom_bar(position=position_fill()) +
	theme_classic() + scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
	labs(x = "Minor Allele Count", y = "Proportion", fill = "Consequence category") + 
	theme(legend.position = "none") +
	scale_fill_brewer(labels = c(
		"Damaging missense",
		"Non-coding",
		"Other missense",
		"PTV",
		"Synonymous",
		"Unknown"), palette="Set3")
print(p)

