library(data.table)
library(dplyr)

# Roel Ophoff data.

df_ucla <- fread("Ophoff_data/table_for_duncan_BIP.txt", colClasses="character")
main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Roel Ophoff")

df_ucla <- df_ucla %>% mutate(SAMPLE_ALIAS = Bipolar_ID)
df_ucla_to_merge <- merge(df_ucla, main_phenotypes, by="SAMPLE_ALIAS") %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, Diagnosis) %>% mutate(PHENOTYPE_FINE_NEW = Diagnosis) %>% select(-Diagnosis)

# Clean-up the naming of the new fine phenotyping

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("Bipolar disorder NOS", "Bipolar NOS")] <- "Bipolar Disorder NOS"

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("Bipolar I")] <- "Bipolar Disorder 1"

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("Bipolar II")] <- "Bipolar Disorder 2"

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("unclear", "Unclear", "unclear, no BP", "exclusion")] <- "Unknown"

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("depression", "Induced BPI", "induced mood disorder", "MDD", "Psychotic disorder NOS")] <- "Other"

df_ucla_to_merge$PHENOTYPE_FINE_NEW[df_ucla_to_merge$PHENOTYPE_FINE_NEW %in%
	c("Schizoaffective disorder")] <- "Schizoaffective"

# And set anything that's missing to it's version in the previous definition

where_missing <- which(df_ucla_to_merge$PHENOTYPE_FINE_NEW == "")
df_ucla_to_merge$PHENOTYPE_FINE_NEW[where_missing] <- df_ucla_to_merge$PHENOTYPE_FINE[where_missing]

# Manually check the differences.
df_ucla_to_merge[df_ucla_to_merge$PHENOTYPE_FINE != df_ucla_to_merge$PHENOTYPE_FINE_NEW,]
# Looks good!
df_ucla_to_merge[is.na(df_ucla_to_merge$PHENOTYPE_FINE_NEW), c("PHENOTYPE_FINE", "PHENOTYPE_FINE_NEW")]

# The new coarse phenotypes for this dataset are the same as the old.
df_ucla_to_merge$PHENOTYPE_COARSE_NEW <- df_ucla_to_merge$PHENOTYPE_COARSE

# Check the column names
names(df_ucla_to_merge)
df_ucla_to_merge <- df_ucla_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# Andrew McQuillin data.

df_ucl <- fread("McQuillin_data/McQuillin_BP_type_data_29Apr2019.txt", colClasses="character")
main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Andrew McQuillin")

names(df_ucl)[which(names(df_ucl) == "Collaborator Sample ID")] <- "SAMPLE_ALIAS" 
names(df_ucl)[which(names(df_ucl) == "Primary Disease")] <- "PHENOTYPE_COARSE_NEW" 

df_ucl_to_merge <- merge(df_ucl, main_phenotypes, by="SAMPLE_ALIAS") %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, BP_type, PHENOTYPE_COARSE_NEW) %>% mutate(PHENOTYPE_FINE_NEW = BP_type) %>% select(-BP_type)

df_ucl_to_merge$PHENOTYPE_FINE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("BP9")] <- "Bipolar Disorder"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("BP1")] <- "Bipolar Disorder 1"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("BP2")] <- "Bipolar Disorder 2"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("SABP")] <- "Schizoaffective"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("SADEP")] <- "Other"
df_ucl_to_merge$PHENOTYPE_COARSE_NEW[(df_ucl_to_merge$PHENOTYPE_COARSE_NEW %in% c("Control", "Controls"))] <- "Control"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[(df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("")) & (df_ucl_to_merge$PHENOTYPE_COARSE_NEW %in% c("Control"))] <- "Control"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[(df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("")) & (df_ucl_to_merge$PHENOTYPE_COARSE_NEW == "Schizophrenia")] <- "Schizophrenia"
df_ucl_to_merge$PHENOTYPE_COARSE_NEW[df_ucl_to_merge$PHENOTYPE_COARSE_NEW %in%
	c("Schizoaffective Disorder", "Schizoaffective Schizophrenia")] <- "Schizoaffective"
df_ucl_to_merge$PHENOTYPE_FINE_NEW[(df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("")) & (df_ucl_to_merge$PHENOTYPE_COARSE_NEW == "Schizoaffective")] <- "Schizoaffective"
df_ucl_to_merge$PHENOTYPE_COARSE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW %in%
	c("Schizoaffective")] <- "Schizoaffective"
df_ucl_to_merge$PHENOTYPE_COARSE_NEW[df_ucl_to_merge$PHENOTYPE_COARSE_NEW %in%
	c("Bipolar Affective Disorder", "Bipolar disorder", "Other PsychosisManic Depressive Disorder",
		"OtherBipolar Psychoses", "OtherMania/Hypomania", "OtherHypomania")] <- "Bipolar Disorder"
df_ucl_to_merge$PHENOTYPE_COARSE_NEW[df_ucl_to_merge$PHENOTYPE_FINE_NEW == "Schizoaffective"] <- "Schizoaffective"

df_ucl_to_merge[df_ucl_to_merge$PHENOTYPE_COARSE != df_ucl_to_merge$PHENOTYPE_COARSE_NEW, c("PHENOTYPE_COARSE", "PHENOTYPE_COARSE_NEW")]
df_ucl_to_merge[df_ucl_to_merge$PHENOTYPE_FINE != df_ucl_to_merge$PHENOTYPE_FINE_NEW, c("PHENOTYPE_FINE", "PHENOTYPE_FINE_NEW")]

# Looks good!
df_ucl_to_merge[is.na(df_ucl_to_merge$PHENOTYPE_FINE_NEW), c("PHENOTYPE_FINE", "PHENOTYPE_FINE_NEW")]

# Check the column names
names(df_ucl_to_merge)
df_ucl_to_merge <- df_ucl_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# Jordan Smoller, Mikael Landen and Nick Craddock data.

# "UK, Cardiff Michael Owen" : icuk...they're actually the Nick Craddock samples.
df_car <- fread("Smoller_data/icuk.output.pheno", colClasses="character") %>% select(IID, bp, bp1, bp2, bp_nos, bp_scz) %>% mutate(PARTICIPANT_ID = IID, PHENOTYPE_COARSE_NEW = bp) %>% select(-c(IID, bp))
df_car$PHENOTYPE_COARSE_NEW[df_car$PHENOTYPE_COARSE_NEW == "1"] <- "Control"
df_car$PHENOTYPE_COARSE_NEW[df_car$PHENOTYPE_COARSE_NEW == "2"] <- "Bipolar Disorder"
df_car$PHENOTYPE_FINE_NEW <- NA
df_car$PHENOTYPE_FINE_NEW[which(df_car$PHENOTYPE_COARSE_NEW == "Control")] <- "Control"
df_car$PHENOTYPE_FINE_NEW[which(df_car$bp1 == "1")] <- "Bipolar Disorder 1"
df_car$PHENOTYPE_FINE_NEW[which(df_car$bp2 == "1")] <- "Bipolar Disorder 2"
df_car$PHENOTYPE_FINE_NEW[which(df_car$bp_nos == "1")] <- "Bipolar Disorder NOS"
df_car$PHENOTYPE_FINE_NEW[which(df_car$bp_scz == "1")] <- "Schizoaffective"

main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Nick Craddock")

# Things that don't match a participant ID in our dataset.
df_car_to_merge <- merge(main_phenotypes, df_car, by="PARTICIPANT_ID", all.x=TRUE)
# It turns out that the data for which Chia-Yen doesn't have subphenotype information, we already have!
df_car_to_merge$PHENOTYPE_FINE[!(df_car_to_merge$PARTICIPANT_ID %in% df_car$PARTICIPANT_ID)]

# Set the phenotypes that are missing to the old versions.
df_car_to_merge$PHENOTYPE_FINE_NEW[is.na(df_car_to_merge$PHENOTYPE_FINE_NEW)] <- df_car_to_merge$PHENOTYPE_FINE[is.na(df_car_to_merge$PHENOTYPE_FINE_NEW)]
df_car_to_merge$PHENOTYPE_COARSE_NEW[is.na(df_car_to_merge$PHENOTYPE_COARSE_NEW)] <- df_car_to_merge$PHENOTYPE_COARSE[is.na(df_car_to_merge$PHENOTYPE_COARSE_NEW)]
df_car_to_merge$PHENOTYPE_COARSE_NEW[df_car_to_merge$PHENOTYPE_FINE_NEW == "Schizoaffective"] <- "Schizoaffective"

df_car_to_merge <- df_car_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE, PHENOTYPE_FINE_NEW)

df_car_to_merge[df_car_to_merge$PHENOTYPE_COARSE != df_car_to_merge$PHENOTYPE_COARSE_NEW, c("PHENOTYPE_COARSE", "PHENOTYPE_COARSE_NEW")]
df_car_to_merge[df_car_to_merge$PHENOTYPE_FINE != df_car_to_merge$PHENOTYPE_FINE_NEW, c("PHENOTYPE_FINE", "PHENOTYPE_FINE_NEW")]

# Looks good!
df_car_to_merge[is.na(df_car_to_merge$PHENOTYPE_FINE_NEW), c("PHENOTYPE_FINE", "PHENOTYPE_FINE_NEW")]

names(df_car_to_merge)

df_car_to_merge <- df_car_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# "USA, MGH Jordan Smoller": mghb
df_mgh <- fread("Smoller_data/mghb.output.pheno", colClasses="character") %>% select(FID, bp, bp1, bp2, bp_nos, bp_scz) %>% mutate(PARTICIPANT_ID = FID, PHENOTYPE_COARSE_NEW = bp) %>% select(-c(FID, bp))
df_mgh$PHENOTYPE_COARSE_NEW[df_mgh$PHENOTYPE_COARSE_NEW == "1"] <- "Control"
df_mgh$PHENOTYPE_COARSE_NEW[df_mgh$PHENOTYPE_COARSE_NEW == "2"] <- "Bipolar Disorder"
df_mgh$PHENOTYPE_FINE_NEW <- NA
df_mgh$PHENOTYPE_FINE_NEW[which(df_mgh$PHENOTYPE_COARSE_NEW == "Control")] <- "Control"
df_mgh$PHENOTYPE_FINE_NEW[which(df_mgh$bp1 == "2")] <- "Bipolar Disorder 1"
df_mgh$PHENOTYPE_FINE_NEW[which(df_mgh$bp2 == "2")] <- "Bipolar Disorder 2"
df_mgh$PHENOTYPE_FINE_NEW[which(df_mgh$bp_nos == "2")] <- "Bipolar Disorder NOS"
df_mgh$PHENOTYPE_FINE_NEW[which(df_mgh$bp_scz == "2")] <- "Schizoaffective"

main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Jordan Smoller")
df_mgh$PARTICIPANT_ID <- gsub(".*\\*", "", df_mgh$PARTICIPANT_ID)

df_mgh_to_merge <- merge(main_phenotypes, df_mgh, by="PARTICIPANT_ID", all.x=TRUE)

# The collection of the Jordan Smoller data for which we don't have subtype information
df_mgh_to_merge$PARTICIPANT_ID[!(df_mgh_to_merge$PARTICIPANT_ID %in% df_mgh$PARTICIPANT_ID)]
# There's still 1000 that are missing.

# The subset that are cases for which we don't have subtype information.
df_mgh_to_merge$PARTICIPANT_ID[!(df_mgh_to_merge$PARTICIPANT_ID %in% df_mgh$PARTICIPANT_ID) & df_mgh_to_merge$PHENOTYPE_COARSE != "Control"]
df_mgh_to_merge$PHENOTYPE_COARSE[!(df_mgh_to_merge$PARTICIPANT_ID %in% df_mgh$PARTICIPANT_ID) & df_mgh_to_merge$PHENOTYPE_COARSE != "Control"]
# There's still 800 cases that Chia-Yen doesn't have a match for.
# Of these, how many don't we have subtype information for?
to_find <- df_mgh_to_merge$PARTICIPANT_ID[!(df_mgh_to_merge$PARTICIPANT_ID %in% df_mgh$PARTICIPANT_ID) & df_mgh_to_merge$PHENOTYPE_FINE == "Bipolar Disorder"]
df_mgh_to_merge$PHENOTYPE_FINE[!(df_mgh_to_merge$PARTICIPANT_ID %in% df_mgh$PARTICIPANT_ID) & df_mgh_to_merge$PHENOTYPE_FINE == "Bipolar Disorder"]

# All of them!
# Found these files on the cluster:
usa1s_pchip_updateid  <- fread("Smoller_data/found_on_cluster/usa1s_pchip_updateid.gencall.pheno", header=FALSE)
usa2smoll_pchip_updateid  <- fread("Smoller_data/found_on_cluster/usa2smoll_pchip_updateid.gencall.pheno", header=FALSE)
usa3con_pchip_updateid  <- fread("Smoller_data/found_on_cluster/usa3con_pchip_updateid.gencall.pheno", header=FALSE)
i2b2_pchip_updateid  <- fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.pheno", header=FALSE)
smol4v11_pchip_updateid  <- fread("Smoller_data/found_on_cluster/smol4v11_pchip_updateid.gencall.pheno", header=FALSE)

colnames(usa1s_pchip_updateid) <- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog")
colnames(usa2smoll_pchip_updateid) <- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog")
colnames(usa3con_pchip_updateid) <- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog")
colnames(i2b2_pchip_updateid) <- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog")
colnames(smol4v11_pchip_updateid) <- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog")

combined <- rbind(usa1s_pchip_updateid, usa2smoll_pchip_updateid, usa3con_pchip_updateid, i2b2_pchip_updateid, smol4v11_pchip_updateid)

# These are the samples that are left. I think stop here.
to_find <- to_find[!(to_find %in% combined$PTID)]

# Check if anything was missed in creating these phenotype files
mgh_txt <- fread("Smoller_data/found_on_cluster/ICCBD.Phenotype.MGH.02Apr14_edit.v2.txt") %>% mutate(IID = as.character(ParticipantID))
usa1s_pchip_updateid_fam  <- fread("Smoller_data/found_on_cluster/usa1s_pchip_updateid.gencall.fam", header=FALSE)
usa2smoll_pchip_updateid_fam  <- fread("Smoller_data/found_on_cluster/usa2smoll_pchip_updateid.gencall.fam", header=FALSE)
usa3con_pchip_updateid_fam  <- fread("Smoller_data/found_on_cluster/usa3con_pchip_updateid.gencall.fam", header=FALSE)
i2b2_pchip_updateid_fam  <- fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.fam", header=FALSE)
smol4v11_pchip_updateid_fam  <- fread("Smoller_data/found_on_cluster/smol4v11_pchip_updateid.gencall.fam", header=FALSE)
mgh_fam <- rbind(usa1s_pchip_updateid_fam, usa2smoll_pchip_updateid_fam,
				 usa3con_pchip_updateid_fam, i2b2_pchip_updateid_fam, smol4v11_pchip_updateid_fam)
colnames(mgh_fam) <- c("PARTICIPANT_ID","IID","PID","MID","gender_fam","pheno_fam")
mgh_check <- merge(mgh_fam, mgh_txt, by="IID")
# None of them are in there. STOP!

# Now that I've found everything I can, need to merge them in.
names(combined)[which(names(combined) == "PTID")] <- "PARTICIPANT_ID"
df_mgh_to_merge <- merge(df_mgh_to_merge, combined, by="PARTICIPANT_ID", all.x=TRUE)
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK <- NA
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_mgh_to_merge$BPTYPE_final == '0')] <- "Control"
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_mgh_to_merge$BPTYPE_final == '1')] <- "Bipolar Disorder 1"
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_mgh_to_merge$BPTYPE_final == '2')] <- "Bipolar Disorder 2"
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_mgh_to_merge$BPTYPE_final == '3')] <- "Bipolar Disorder NOS"
df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_mgh_to_merge$BPTYPE_final == '4')] <- "Schizoaffective"

# Ensure that PHENOTYPE_FINE_NEW_CHECK and PHENOTYPE_FINE_NEW match where they are defined.
sum(df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK != df_mgh_to_merge$PHENOTYPE_FINE_NEW, na.rm=TRUE)
# Good.

# Now set those phenotypes that are missing in PHENOTYPE_FINE_NEW to the non-missing values in PHENOTYPE_FINE_NEW_CHECK.
where <- which(is.na(df_mgh_to_merge$PHENOTYPE_FINE_NEW) & !is.na(df_mgh_to_merge$PHENOTYPE_FINE_NEW_CHECK))
df_mgh_to_merge[where, "PHENOTYPE_FINE_NEW"] <- df_mgh_to_merge[where, "PHENOTYPE_FINE_NEW_CHECK"]
# Everything that is missing in this final fine grained phenotyping should now get set to the previous 'FINE' phenotype (which is just coarse).

where <- which(is.na(df_mgh_to_merge$PHENOTYPE_FINE_NEW))
df_mgh_to_merge[where, "PHENOTYPE_FINE_NEW"] <- df_mgh_to_merge[where, "PHENOTYPE_FINE"]

# Set the coarse and fine new values that are undefined to the old values
df_mgh_to_merge$PHENOTYPE_COARSE_NEW[is.na(df_mgh_to_merge$PHENOTYPE_COARSE_NEW)] <- df_mgh_to_merge$PHENOTYPE_COARSE[is.na(df_mgh_to_merge$PHENOTYPE_COARSE_NEW)]
df_mgh_to_merge$PHENOTYPE_COARSE_NEW[df_mgh_to_merge$PHENOTYPE_FINE_NEW == "Schizoaffective"] <- "Schizoaffective"

# Remove the duplicated row
df_mgh_to_merge <- df_mgh_to_merge[-duplicated(df_mgh_to_merge$SAMPLE_ALIAS),]
df_mgh_to_merge <- df_mgh_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# "Karolinska Mikael Landen" : swei, swa2, lah1
# df_kl1 <- fread("Smoller_data/swa2.output.pheno", colClasses="character") %>% mutate(PARTICIPANT_ID = FID, PHENOTYPE_COARSE_NEW = bp) %>% select(-c(FID, bp))
# df_kl1_fam <- fread("Smoller_data/bip_swa2_eur_sr-qc.fam")
# colnames(df_kl1_fam) <- c("FID","PARTICIPANT_ID","PID","MID","gender_fam","pheno_fam")
# df_kl1_covar <- fread("Smoller_data/swa2.output.covar") %>% mutate(PARTICIPANT_ID = FID)
# df_kl1 <- merge(df_kl1, df_kl1_covar, by="PARTICIPANT_ID")
# df_kl1$PARTICIPANT_ID <- gsub(".*\\*", "", df_kl1$PARTICIPANT_ID)
# df_kl1 <- merge(df_kl1, df_kl1_fam, by="PARTICIPANT_ID", all=TRUE)

# df_kl2 <- fread("Smoller_data/swei.output.pheno", colClasses="character") %>% mutate(PARTICIPANT_ID = FID, PHENOTYPE_COARSE_NEW = bp) %>% select(-c(FID, bp))
# df_kl2_fam <- fread("Smoller_data/bip_swei_eur_sr-qc.fam")
# colnames(df_kl2_fam) <- c("FID","PARTICIPANT_ID","PID","MID","gender_fam","pheno_fam")
# df_kl2_covar <- fread("Smoller_data/swei.output.covar") %>% mutate(PARTICIPANT_ID = FID)
# df_kl2 <- merge(df_kl2, df_kl2_covar, by="PARTICIPANT_ID")
# df_kl2$PARTICIPANT_ID <- gsub(".*\\*", "", df_kl2$PARTICIPANT_ID)
# df_kl2 <- merge(df_kl2, df_kl2_fam, by="PARTICIPANT_ID", all=TRUE)

# df_kl3 <- fread("Smoller_data/lah1.output.pheno", colClasses="character") %>% mutate(PARTICIPANT_ID = FID, PHENOTYPE_COARSE_NEW = bp) %>% select(-c(FID, bp))
# df_kl3_fam <- fread("Smoller_data/bip_lah1_eur_rk-qc.fam")
# colnames(df_kl3_fam) <- c("PARTICIPANT_ID","FID","PID","MID","gender_fam","pheno_fam")
# df_kl3_covar <- fread("Smoller_data/lah1.output.covar") %>% mutate(PARTICIPANT_ID = FID)
# df_kl3 <- merge(df_kl3, df_kl3_covar, by="PARTICIPANT_ID")
# df_kl3$PARTICIPANT_ID <- gsub(".*\\*", "", df_kl3$PARTICIPANT_ID)
# df_kl3_fam$PARTICIPANT_ID <- gsub(".*\\*", "", df_kl3_fam$PARTICIPANT_ID)
# df_kl3 <- merge(df_kl3, df_kl3_fam, by="PARTICIPANT_ID", all=TRUE)

# Latest information from Karolinska.
df_kl <- fread("Landen_data/swedish_bipolar_subtype.csv", sep=',')
df_kl$PHENOTYPE_FINE_NEW <- NA
df_kl$PHENOTYPE_FINE_NEW[which(df_kl$BP1 == 1)] <- "Bipolar Disorder 1"
df_kl$PHENOTYPE_FINE_NEW[which(df_kl$BP2 == 1)] <- "Bipolar Disorder 2"
df_kl$PHENOTYPE_FINE_NEW[which(df_kl$BP_NOS == 1)] <- "Bipolar Disorder NOS"
df_kl$PHENOTYPE_FINE_NEW[which(df_kl$SCZ_BP == 1)] <- "Schizoaffective"
df_kl$PHENOTYPE_FINE_NEW[which(df_kl$ROLE == "Control")] <- "Control"
df_kl$PHENOTYPE_FINE_NEW[(df_kl$ROLE == 'Case') & (is.na(df_kl$PHENOTYPE_FINE_NEW))] <- "Bipolar Disorder"

df_kl$PHENOTYPE_COARSE_NEW <- NA
df_kl$PHENOTYPE_COARSE_NEW[which(df_kl$PHENOTYPE_FINE_NEW == "Control")] <- "Control"
df_kl$PHENOTYPE_COARSE_NEW[which(df_kl$PHENOTYPE_FINE_NEW %in% c("Bipolar Disorder 1", "Bipolar Disorder 2", "Bipolar Disorder", "Bipolar Disorder NOS"))] <- "Bipolar Disorder" 
df_kl$PHENOTYPE_COARSE_NEW[which(df_kl$PHENOTYPE_FINE_NEW == "Schizoaffective")] <- "Schizoaffective"

main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Mikael Landen")
df_kl_to_merge <- merge(df_kl, main_phenotypes, by='PARTICIPANT_ID', all=TRUE)

# Check to ensure that all the FIDs are unique
# any(duplicated(c(df_kl1$PARTICIPANT_ID, df_kl2$PARTICIPANT_ID, df_kl3$PARTICIPANT_ID)))
# Nope, good.

# df_kl <- rbind(df_kl1, df_kl2, df_kl3)

# df_kl$PHENOTYPE_COARSE_NEW[df_kl$PHENOTYPE_COARSE_NEW == "1"] <- "Control"
# df_kl$PHENOTYPE_COARSE_NEW[df_kl$PHENOTYPE_COARSE_NEW == "2"] <- "Bipolar Disorder"
# df_kl$PHENOTYPE_FINE_NEW <- NA
# df_kl$PHENOTYPE_FINE_NEW[which(df_kl$PHENOTYPE_COARSE_NEW == "Control")] <- "Control"
# df_kl$PHENOTYPE_FINE_NEW[which(df_kl$bp1 == "1")] <- "Bipolar Disorder 1"
# df_kl$PHENOTYPE_FINE_NEW[which(df_kl$bp2 == "1")] <- "Bipolar Disorder 2"
# df_kl$PHENOTYPE_FINE_NEW[which(df_kl$bp_nos == "1")] <- "Bipolar Disorder"
# df_kl$PHENOTYPE_FINE_NEW[which(df_kl$bp_scz == "1")] <- "Schizoaffective"

# df_kl_to_merge <- merge(main_phenotypes, df_kl, by="PARTICIPANT_ID", all.x=TRUE)
# df_kl_to_merge$gender[df_kl_to_merge$gender == "1"] <- "Male"
# df_kl_to_merge$gender[df_kl_to_merge$gender == "2"] <- "Female"

# df_kl_to_merge$pheno_fam[df_kl_to_merge$pheno_fam == 1] <- "Control"
# df_kl_to_merge$pheno_fam[df_kl_to_merge$pheno_fam == 2] <- "Bipolar Disorder"
# df_kl_to_merge$pheno_fam[df_kl_to_merge$pheno_fam == -9] <- NA

# # There's a lot of differences at the subtype level.
# dim(df_kl_to_merge[which(df_kl_to_merge$PHENOTYPE_FINE == df_kl_to_merge$PHENOTYPE_FINE_NEW),])
# dim(df_kl_to_merge[which(df_kl_to_merge$PHENOTYPE_FINE != df_kl_to_merge$PHENOTYPE_FINE_NEW),])

# # Are there the same amount of differences at the case level?
# dim(df_kl_to_merge[which(df_kl_to_merge$PHENOTYPE_COARSE != df_kl_to_merge$PHENOTYPE_COARSE_NEW),])
# df_kl_to_merge[which(df_kl_to_merge$PHENOTYPE_COARSE != df_kl_to_merge$PHENOTYPE_COARSE_NEW), c("SAMPLE_ALIAS", "PARTICIPANT_ID", "PHENOTYPE_COARSE", "PHENOTYPE_COARSE_NEW")]

# df_kl_to_merge$PARTICIPANT_ID[!(df_kl_to_merge$PARTICIPANT_ID %in% df_kl$PARTICIPANT_ID)]

# # There's still 1000 that are missing.

# # The subset that are cases for which we don't have subtype information.
# df_kl_to_merge$PARTICIPANT_ID[!(df_kl_to_merge$PARTICIPANT_ID %in% df_kl$PARTICIPANT_ID) & df_kl_to_merge$PHENOTYPE_COARSE != "Control"]
# df_kl_to_merge$PHENOTYPE_COARSE[!(df_kl_to_merge$PARTICIPANT_ID %in% df_kl$PARTICIPANT_ID) & df_kl_to_merge$PHENOTYPE_COARSE != "Control"]
# # There's still 400 cases that Chia-Yen doesn't have a match for.
# # Of these, how many don't we have subtype information for?
# to_find <- df_kl_to_merge$PARTICIPANT_ID[!(df_kl_to_merge$PARTICIPANT_ID %in% df_kl$PARTICIPANT_ID) & df_kl_to_merge$PHENOTYPE_FINE == "Bipolar Disorder"]
# df_kl_to_merge$PHENOTYPE_FINE[!(df_kl_to_merge$PARTICIPANT_ID %in% df_kl$PARTICIPANT_ID) & df_kl_to_merge$PHENOTYPE_FINE == "Bipolar Disorder"]
# # 143 are left. Let's look for these.

# # Let's read in the Swedish information that I have and see if it's in there.
# swe_cluster_df <- fread("Smoller_data/found_on_cluster/report_phenotype_sweden_07JAN2015_edit.txt")
# car_cluster_df <- fread("Smoller_data/found_on_cluster/ICCBD.Phenotype.Cardiff.October.07.2014.txt")

# to_find <- to_find[!(to_find %in% swe_cluster_df$participant_id)]

# names(swe_cluster_df)[which(names(swe_cluster_df) == "participant_id")] <- "PARTICIPANT_ID"
# df_kl_to_merge <- merge(df_kl_to_merge, swe_cluster_df, by="PARTICIPANT_ID", all.x=TRUE)
# df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK <- NA
# df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_kl_to_merge$BP1 == '1')] <- "Bipolar Disorder 1"
# df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_kl_to_merge$BP2 == '1')] <- "Bipolar Disorder 2"
# df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_kl_to_merge$BP_NOS == '1')] <- "Bipolar Disorder"
# df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK[which(df_kl_to_merge$SCZ_BP == '1')] <- "Schizoaffective"

# df_kl_to_merge$PHENOTYPE_COARSE_NEW_CHECK <- df_kl_to_merge$pheno_fam
# any(df_kl_to_merge$PHENOTYPE_COARSE_NEW != df_kl_to_merge$PHENOTYPE_COARSE_NEW_CHECK)
# # None - good.

# any(df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK != df_kl_to_merge$PHENOTYPE_FINE_NEW)
# # None - good

# where <- which(!is.na(df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK))
# df_kl_to_merge$PHENOTYPE_FINE_NEW[where] <- df_kl_to_merge$PHENOTYPE_FINE_NEW_CHECK[where]

# where <- which(!is.na(df_kl_to_merge$PHENOTYPE_COARSE_NEW_CHECK))
# df_kl_to_merge$PHENOTYPE_COARSE_NEW[where] <- df_kl_to_merge$PHENOTYPE_COARSE_NEW_CHECK[where]

# # 2 left! Stop here.
# # Plan: Set the Bipolar subtype to Chia-Yen's bipolar subtype if it is available.
# # If it is not available, set the subtype to the current PHENOTYPE_FINE data.
# where <- which(is.na(df_kl_to_merge$PHENOTYPE_COARSE_NEW))
# df_kl_to_merge$PHENOTYPE_COARSE_NEW[where] <- df_kl_to_merge$PHENOTYPE_COARSE[where]

# where <- which(is.na(df_kl_to_merge$PHENOTYPE_FINE_NEW))
# df_kl_to_merge$PHENOTYPE_FINE_NEW[where] <- df_kl_to_merge$PHENOTYPE_FINE[where]
# df_kl_to_merge <- df_kl_to_merge[-which(duplicated(df_kl_to_merge$PARTICIPANT_ID)),]

df_kl_to_merge <- df_kl_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)
# Andreas Reif data
df_ger <- fread("Reif_data/Reif_BD_Info_toDuncan_672019_TMK.txt")
main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Andreas Reif")

df_ger <- merge(df_ger, main_phenotypes, by='SAMPLE_ALIAS', all=TRUE)
names(df_ger)[names(df_ger) == "Bipolar II (I/II/NOS/U)"] <- "PHENOTYPE_FINE_NEW"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_COARSE == "Control")] <- "Control"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_FINE_NEW == "Bipolar-III")] <- "Bipolar Disorder"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_FINE_NEW == "Bipolar-II")] <- "Bipolar Disorder 2"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_FINE_NEW == "Bipolar-I")] <- "Bipolar Disorder 1"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_FINE_NEW == "Bipolar NOS")] <- "Bipolar Disorder NOS"
df_ger$PHENOTYPE_FINE_NEW[which(df_ger$PHENOTYPE_FINE_NEW == "U")] <- "Unknown"
df_ger$PHENOTYPE_FINE_NEW[is.na(df_ger$PHENOTYPE_FINE_NEW)] <- "Unknown"

# Send unknowns in the subtypes to unknown in the coarse phenotype.
df_ger$PHENOTYPE_COARSE_NEW <- df_ger$PHENOTYPE_FINE_NEW
df_ger$PHENOTYPE_COARSE_NEW[grep("Bipolar", df_ger$PHENOTYPE_FINE_NEW)] <- "Bipolar Disorder"

df_ger_to_merge <- df_ger %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# Bob Yolken data

phenotypes <- fread("Yolken_data/seph_phenotyping_05_20_19_dp.txt", colClasses="character")
DSMIV <- fread("Yolken_data/DSM-IV.txt", colClasses="character")
names(DSMIV) <- c("dsmivdiag", "DSM_translation")

merged_df <- merge(phenotypes, DSMIV, by = "dsmivdiag", all.x=TRUE)

merged_df %>% filter(is.na(DSM_translation))

merged_df$DSM_translation <- gsub("\xca", "", merged_df$DSM_translation)
merged_df$primary_disease_subtype <- merged_df$DSM_translation
merged_df$primary_disease_subtype[grep("Bipolar I Disorder", merged_df$primary_disease_subtype)] <- "Bipolar Disorder 1"
merged_df$primary_disease_subtype[grep("Bipolar II Disorder", merged_df$primary_disease_subtype)] <- "Bipolar Disorder 2"
merged_df$primary_disease_subtype[grep("Bipolar Disorder NOS", merged_df$primary_disease_subtype)] <- "Bipolar Disorder NOS"
merged_df$primary_disease_subtype[grep("Major Depressive Disorder", merged_df$primary_disease_subtype)] <- "Other"

merged_df <- merged_df %>% mutate(SAMPLE_ALIAS = ry_id) %>% mutate(PHENOTYPE_FINE_NEW = primary_disease_subtype) %>% select(SAMPLE_ALIAS, PHENOTYPE_FINE_NEW)
merged_df$PHENOTYPE_FINE_NEW[which(is.na(merged_df$PHENOTYPE_FINE_NEW))] <- "Control"

# Make sure that this matches up with the names in the main phenotype file.
main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv") %>% filter(PI == "Bob Yolken")
df_jh_to_merge <- merge(main_phenotypes, merged_df, by="SAMPLE_ALIAS", all.x=TRUE)
df_jh_to_merge$PHENOTYPE_FINE_NEW[is.na(df_jh_to_merge$PHENOTYPE_FINE_NEW)] <- df_jh_to_merge$PHENOTYPE_FINE[is.na(df_jh_to_merge$PHENOTYPE_FINE_NEW)]

df_jh_to_merge$PHENOTYPE_COARSE_NEW <- df_jh_to_merge$PHENOTYPE_FINE_NEW
df_jh_to_merge$PHENOTYPE_COARSE_NEW[grep("Bipolar", df_jh_to_merge$PHENOTYPE_COARSE_NEW)] <- "Bipolar Disorder"

df_jh_to_merge <- df_jh_to_merge %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE, PHENOTYPE_FINE, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# Now, combine everything

df_merged <- rbind(df_ucla_to_merge, df_ucl_to_merge, df_car_to_merge, df_mgh_to_merge, df_kl_to_merge, df_ger_to_merge, df_jh_to_merge)
# Remove the few stray duplicates

df_merged <- df_merged[-which(duplicated(df_merged$SAMPLE_ALIAS)),] %>% select(SAMPLE_ALIAS, PHENOTYPE_COARSE_NEW, PHENOTYPE_FINE_NEW)

# Now, merge this with the original data, and set the new coarse and fine information to the old where it is missing (primarily Edinburgh).
main_phenotypes <- fread("BIP_phenotype_information_cleaned.tsv")
main_phenotypes_merged <- merge(main_phenotypes, df_merged, by='SAMPLE_ALIAS', all=TRUE)

where_coarse_new_missing <- which(is.na(main_phenotypes_merged$PHENOTYPE_COARSE_NEW))
main_phenotypes_merged$PHENOTYPE_COARSE_NEW[where_coarse_new_missing] <- main_phenotypes_merged$PHENOTYPE_COARSE[where_coarse_new_missing]

where_fine_new_missing <- which(is.na(main_phenotypes_merged$PHENOTYPE_FINE_NEW))
main_phenotypes_merged$PHENOTYPE_FINE_NEW[where_fine_new_missing] <- main_phenotypes_merged$PHENOTYPE_FINE[where_fine_new_missing]

# Remove first row - it's just an empty line...
main_phenotypes_merged <- main_phenotypes_merged[-1,]

fwrite(main_phenotypes_merged, "BIP_phenotype_information_cleaned_new_subtype_information_added.tsv", sep='\t')
