library(data.table)
library(dplyr)

dt_pheno <- fread("../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added.tsv")

dt_pheno$PHENOTYPE_FINE_NEW_CHECK <- dt_pheno$PHENOTYPE_FINE_NEW
dt_pheno$PHENOTYPE_FINE_NEW_CHECK[grep("Bipolar Disorder", dt_pheno$PHENOTYPE_FINE_NEW_CHECK)] <- "Bipolar Disorder"
# These are the mismatches:
# dt_pheno[which(dt_pheno$PHENOTYPE_FINE_NEW_CHECK != dt_pheno$PHENOTYPE_COARSE_NEW),]
dt_pheno$PHENOTYPE_FINE_NEW[(dt_pheno$PHENOTYPE_FINE_NEW_CHECK == "Unknown") & (dt_pheno$PHENOTYPE_COARSE_NEW == "Other")] <- "Other"
dt_pheno$PHENOTYPE_FINE_NEW_CHECK[(dt_pheno$PHENOTYPE_FINE_NEW_CHECK == "Unknown") & (dt_pheno$PHENOTYPE_COARSE_NEW == "Other")] <- "Other"
dt_pheno$PHENOTYPE_FINE_NEW[which(dt_pheno$PHENOTYPE_FINE_NEW_CHECK != dt_pheno$PHENOTYPE_COARSE_NEW)] <- "Bipolar Disorder"
fwrite(dt_pheno, file="../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised.tsv")
