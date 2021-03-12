library(data.table)
library(dplyr)

# Counts for supplementary tables.
dt <- fread("~/Repositories/BipEx/phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv")

# Table S3, before filtering to QCed samples.
print(dt %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_FINE, LOCATION) %>% summarise(sum(PSYCHOSIS), sum(!PSYCHOSIS)), n=Inf)
# Totals
# By fine-grain phenotype
print(dt %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_FINE) %>% summarise(sum(PSYCHOSIS),  sum(!PSYCHOSIS)))
# By Location
print(dt %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_COARSE, LOCATION) %>% count())
print(dt %>% filter(!is.na(PSYCHOSIS), PHENOTYPE_COARSE == "Bipolar Disorder") %>% count())

# Table S1, counts before QC by PI.
print(dt %>% group_by(PHENOTYPE_FINE, PI, LOCATION) %>% summarise(n()), n=Inf)
print(dt %>% group_by(PHENOTYPE_FINE) %>% summarise(n()), n=Inf)


# Following restriction to high quality sequence data.
# This is the output on the cloud, following all QC.

# Table S4
SAMPLE_AFTER_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz | gzcat' 
dt_qc <- fread(cmd = SAMPLE_AFTER_QC_FILE)
print(dt_qc %>% filter(!(PHENOTYPE_FINE %in% c("Schizophrenia", "Unknown", "Other"))) %>% group_by(PHENOTYPE_FINE, LOCATION) %>% count(), n=Inf)
print(dt_qc %>% filter(!(PHENOTYPE_COARSE %in% c("Schizophrenia", "Unknown", "Other"))) %>% group_by(PHENOTYPE_COARSE, LOCATION) %>% count(), n=Inf)

# Totals
print(dt_qc %>% filter(!(PHENOTYPE_FINE %in% c("Schizophrenia", "Unknown", "Other"))) %>% group_by(PHENOTYPE_FINE) %>% count(), n=Inf)
print(dt_qc %>% filter(!(PHENOTYPE_COARSE %in% c("Schizophrenia", "Unknown", "Other"))) %>% group_by(PHENOTYPE_COARSE) %>% count(), n=Inf)

print(dt_qc %>% filter(PHENOTYPE_COARSE %in% c("Bipolar Disorder", "Control")) %>% group_by(LOCATION) %>% count())
print(dt_qc %>% filter(PHENOTYPE_COARSE %in% c("Bipolar Disorder", "Control")) %>% count())

# Table S3, after filtering to QCed samples.
print(dt_qc %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_FINE, LOCATION) %>% summarise(sum(PSYCHOSIS), sum(!PSYCHOSIS)), n=Inf)
# Totals
# By fine-grain phenotype
print(dt_qc %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_FINE) %>% summarise(sum(PSYCHOSIS),  sum(!PSYCHOSIS)))
# By Location
print(dt_qc %>% filter(!is.na(PSYCHOSIS)) %>% group_by(PHENOTYPE_COARSE, LOCATION) %>% count())
print(dt_qc %>% filter(!is.na(PSYCHOSIS), PHENOTYPE_COARSE == "Bipolar Disorder") %>% count())

# Create small barplots for Figure 1
dt_qc <- fread(cmd = SAMPLE_AFTER_QC_FILE)
dt_bar <- as.data.table(dt_qc %>% filter(!(PHENOTYPE_FINE %in% c("Schizophrenia", "Unknown", "Other", "Schizoaffective"))) %>% group_by(PHENOTYPE_FINE) %>% summarise(count=n()))
dt_bar <- dt_bar %>% mutate(PHENOTYPE_FINE = as.factor(PHENOTYPE_FINE))
dt_bar <- dt_bar %>% mutate(PHENOTYPE_COARSE = ifelse(PHENOTYPE_FINE == "Control", "Control", "Bipolar Disorder"))
pdf("Figure_1_barplots.pdf")
ggplot(dt_bar, aes(fill=PHENOTYPE_FINE, y=count, x=PHENOTYPE_COARSE)) + geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(dt_bar %>% filter(PHENOTYPE_FINE %in% c("Bipolar Disorder 1", "Control")), aes(fill=PHENOTYPE_FINE, y=count, x=PHENOTYPE_COARSE)) + geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(dt_bar %>% filter(PHENOTYPE_FINE %in% c("Bipolar Disorder 2", "Control")), aes(fill=PHENOTYPE_FINE, y=count, x=PHENOTYPE_COARSE)) + geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()