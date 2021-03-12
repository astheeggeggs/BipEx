library(data.table)
library(dplyr)

# First, we now have even more information relating to the phenotypes from the ICCBD projects.
# Let's ensure that the fine and coarse phenotyping matches.

# Cardiff 
cat("Clean and merge in Cardiff information...\n\n")
dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv") %>% filter(LOCATION == "Cardiff, UK")

# The following two files should contain the same data.

# ICCBD data
dt_ICCBD_cardiff <- fread("ICCBD_data/ICCBD.Phenotype.Cardiff.October.07.2014.txt") %>% rename(PARTICIPANT_ID=ID) %>%
	mutate(PHENOTYPE_COARSE = case_when(SCZ_BP == 1 ~ "Schizoaffective", 
			Type == 1 ~ "Bipolar Disorder",
			Type == 0 ~ "Control"),
		PHENOTYPE_FINE = case_when(
			BP1 == 1 ~ "Bipolar Disorder 1",
			BP2 == 1 ~ "Bipolar Disorder 2",
			BP_NOS == 1 ~ "Bipolar Disorder NOS",
			SCZ_BP == 1 ~ "Schizoaffective",
			Type == 1 ~ "Bipolar Disorder",
			TRUE ~ NA_character_)
		) %>% select(PARTICIPANT_ID, PHENOTYPE_COARSE, PHENOTYPE_FINE, PSYCHOSIS = Psychosis,
			# Adding in the age of onset variables as well
			# Age of first impairment, binned.
			AGE_FI_24 = Age_FI_24, AGE_FI_40 = Age_FI_40,
			# Age of first symptoms, binned.
			AGE_FS_24 = Age_FS_24, AGE_FS_40 = Age_FS_40,
			# Age of first diagnosis, binned (Note that this is NA for all the Cardiff data).
			AGE_D_24 = Age_D_24, AGE_D_40 = Age_D_40
		)

# Data grabbed from Chia-yen's folder
dt_smoller_cardiff <- fread("Smoller_data/icuk.output.pheno", colClasses="character") %>% 
select(
	PARTICIPANT_ID = IID, PHENOTYPE_COARSE = bp, PSYCHOSIS=psychosis, 
	bp1, bp2, bp_nos, bp_scz, 
	# Adding in the age of onset variables as well
	# Age of first impairment, binned.
	AGE_FI_24 = age_fi_24, AGE_FI_40 = age_fi_40,
	# Age of first symptoms, binned.
	AGE_FS_24 = age_fs_24, AGE_FS_40 = age_fs_40,
	# Age of first diagnosis, binned (Note that this is NA for all the Cardiff data).
	AGE_D_24 = age_d_24, AGE_D_40 = age_d_40
)

dt_smoller_cardiff <- dt_smoller_cardiff %>% 
	mutate(PHENOTYPE_COARSE = 
		case_when(
			bp_scz == "1" ~ "Schizoaffective",
			PHENOTYPE_COARSE == "2" ~ "Bipolar Disorder", 
			PHENOTYPE_COARSE == "1" ~ "Control",
			TRUE ~ NA_character_)
		) %>% 
	mutate(PHENOTYPE_FINE = 
		case_when(
			PHENOTYPE_COARSE == "Control" ~ "Control",
			bp1 == "1" ~ "Bipolar Disorder 1",
			bp2 == "1" ~ "Bipolar Disorder 2",
			bp_nos == "1" ~ "Bipolar Disorder NOS",
			bp_scz == "1" ~ "Schizoaffective",
			TRUE ~ NA_character_)
		)

# They all match at the intersection. Great.
dim(merge(dt_ICCBD_cardiff, dt_smoller_cardiff, by=c("PARTICIPANT_ID", "PHENOTYPE_COARSE", "PHENOTYPE_FINE"))) == dim(merge(dt_ICCBD_cardiff, dt_smoller_cardiff, by="PARTICIPANT_ID"))[1]

# Now, grab all the Cardiff information from these two merged files.
dt_cardiff <- merge(dt_ICCBD_cardiff, dt_smoller_cardiff, by=c("PARTICIPANT_ID", "PHENOTYPE_COARSE", "PHENOTYPE_FINE"), all=TRUE)
# Send stupid values to NA
dt_cardiff <- dt_cardiff %>% 
	mutate(PSYCHOSIS.x = ifelse(PSYCHOSIS.x %in% c(0,1), PSYCHOSIS.x, NA),
		   PSYCHOSIS.y = ifelse(PSYCHOSIS.y %in% c(0,1), PSYCHOSIS.y, NA),
		   # Age of first impairment, binned.
		   AGE_FI_24.x = ifelse(AGE_FI_24.x %in% c(0,1,2), AGE_FI_24.x, NA),
		   AGE_FI_24.y = ifelse(AGE_FI_24.y %in% c(0,1,2), AGE_FI_24.y, NA),
		   AGE_FI_40.x = ifelse(AGE_FI_40.x %in% c(0,1,2), AGE_FI_40.x, NA),
		   AGE_FI_40.y = ifelse(AGE_FI_40.y %in% c(0,1,2), AGE_FI_40.y, NA),
		   # Age of first symptoms, binned.
		   AGE_FS_24.x = ifelse(AGE_FS_24.x %in% c(0,1,2), AGE_FS_24.x, NA),
		   AGE_FS_24.y = ifelse(AGE_FS_24.y %in% c(0,1,2), AGE_FS_24.y, NA),
		   AGE_FS_40.x = ifelse(AGE_FS_40.x %in% c(0,1,2), AGE_FS_40.x, NA),
		   AGE_FS_40.y = ifelse(AGE_FS_40.y %in% c(0,1,2), AGE_FS_40.y, NA),
		   # Age of first diagnosis, binned (Note that this is NA for all the Cardiff data).
		   AGE_D_24.x = ifelse(AGE_D_24.x %in% c(0,1,2), AGE_D_24.x, NA),
		   AGE_D_24.y = ifelse(AGE_D_24.y %in% c(0,1,2), AGE_D_24.y, NA),
		   AGE_D_40.x = ifelse(AGE_D_40.x %in% c(0,1,2), AGE_D_40.x, NA),
		   AGE_D_40.y = ifelse(AGE_D_40.y %in% c(0,1,2), AGE_D_40.y, NA)
		   )

# Ensure that the psychosis information matches.
if (all(dt_cardiff$PSYCHOSIS.x == dt_cardiff$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

# Ensure that all the age of onset information matches.

# Impairment
if (all(dt_cardiff$AGE_FI_24.x == dt_cardiff$AGE_FI_24.y, na.rm=TRUE)) {
	cat("Age of onset; first impairment (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first impairment (24) does not agree in these two files\n")
}

if (all(dt_cardiff$AGE_FI_40.x == dt_cardiff$AGE_FI_40.y, na.rm=TRUE)) {
	cat("Age of onset; first impairment (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first impairment (40) does not agree in these two files\n")
}

# First symptoms
if (all(dt_cardiff$AGE_FS_24.x == dt_cardiff$AGE_FS_24.y, na.rm=TRUE)) {
	cat("Age of onset; first symptoms (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first symptoms (24) does not agree in these two files\n")
}

if (all(dt_cardiff$AGE_FS_40.x == dt_cardiff$AGE_FS_40.y, na.rm=TRUE)) {
	cat("Age of onset; first symptoms (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first symptoms (40) does not agree in these two files\n")
}

# Diagnosis
if (all(dt_cardiff$AGE_D_24.x == dt_cardiff$AGE_D_24.y, na.rm=TRUE)) {
	cat("Age of onset; diagnosis (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; diagnosis (24) does not agree in these two files\n")
}

if (all(dt_cardiff$AGE_D_40.x == dt_cardiff$AGE_D_40.y, na.rm=TRUE)) {
	cat("Age of onset; diagnosis (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; diagnosis (40) does not agree in these two files\n")
}

# Where one of the values is missing, pick the non-missing value
dt_cardiff <- dt_cardiff %>% 
	mutate(
		PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x),
		AGE_FI_24 = ifelse(is.na(AGE_FI_24.x), AGE_FI_24.y, AGE_FI_24.x),
		AGE_FI_40 = ifelse(is.na(AGE_FI_40.x), AGE_FI_40.y, AGE_FI_40.x),
		AGE_FS_24 = ifelse(is.na(AGE_FS_24.x), AGE_FS_24.y, AGE_FS_24.x),
		AGE_FS_40 = ifelse(is.na(AGE_FS_40.x), AGE_FS_40.y, AGE_FS_40.x),
		AGE_D_24 = ifelse(is.na(AGE_D_24.x), AGE_D_24.y, AGE_D_24.x),
		AGE_D_40 = ifelse(is.na(AGE_D_40.x), AGE_D_40.y, AGE_D_40.x)
		)

# These also all match at the intersection. Great.
dt <- merge(dt,
	(dt_cardiff %>% 
		mutate(
			AGE_FI_24 = as.numeric(AGE_FI_24),
			AGE_FI_40 = as.numeric(AGE_FI_40),
			AGE_FS_24 = as.numeric(AGE_FS_24),
			AGE_FS_40 = as.numeric(AGE_FS_40),
			AGE_D_24 = as.numeric(AGE_D_24),
			AGE_D_40 = as.numeric(AGE_D_40),
			PSYCHOSIS = as.numeric(PSYCHOSIS)) %>% 
		select(
			PARTICIPANT_ID, PHENOTYPE_COARSE, PHENOTYPE_FINE, PSYCHOSIS,
			# Adding in the age of onset variables as well
			AGE_FI_24, AGE_FI_40, AGE_FS_24, AGE_FS_40, AGE_D_24, AGE_D_40
		)
	), by=c("PARTICIPANT_ID", "PHENOTYPE_COARSE", "PHENOTYPE_FINE"), all.x=TRUE)

# First, double check that the data that Di-Florio sent over matches in terms of Bipolar case-control information

# Data sent in February 2020
dt_di_florio_cardiff <- fread("Di_Florio_data/Cardiff_Phenotype_Query_17022020.tsv", header=TRUE)
dt_di_florio_cardiff <- dt_di_florio_cardiff %>% mutate(
	PHENOTYPE_FINE = 
		ifelse(Diagnosis == "BPI", "Bipolar Disorder 1",
			ifelse(Diagnosis == "BPII", "Bipolar Disorder 2",
				ifelse(Diagnosis == "BP NOS", "Bipolar Disorder NOS",
					ifelse(Diagnosis == "SA BP", "Schizoaffective", NA)
				)

			)
		)
	)

dt_di_florio_cardiff <- dt_di_florio_cardiff %>% mutate(PHENOTYPE_COARSE = ifelse(grepl("BP", Diagnosis), "Bipolar Disorder", NA),
	PARTICIPANT_ID = V1, PSYCHOSIS = Psychosis) %>% select(PARTICIPANT_ID, PHENOTYPE_COARSE, PHENOTYPE_FINE, PSYCHOSIS)

dt_di_florio_cardiff <- merge(dt, dt_di_florio_cardiff, by = "PARTICIPANT_ID", all.y=TRUE)

# Check agreement of phenotypes and subtypes where they are defined.
if (all(dt_di_florio_cardiff$PHENOTYPE_COARSE.x == dt_di_florio_cardiff$PHENOTYPE_COARSE.y, na.rm=TRUE)) {
	cat("Coarse phenotype information (where defined) agrees between these two files\n")
} else {
	cat("Coarse phenotype information does not agree in these two files\n")
}

if (all(dt_di_florio_cardiff$PHENOTYPE_FINE.x == dt_di_florio_cardiff$PHENOTYPE_FINE.y, na.rm=TRUE)) {
	cat("Fine phenotype information (where defined) agrees between these two files\n")
} else {
	cat("Fine phenotype information does not agree in these two files\n")
	print(dt_di_florio_cardiff %>% filter(dt_di_florio_cardiff$PHENOTYPE_FINE.x != dt_di_florio_cardiff$PHENOTYPE_FINE.y) %>% select(PHENOTYPE_FINE.x, PHENOTYPE_FINE.y))
	# Change the result to what's been recently sent.
	samples_different <- dt_di_florio_cardiff %>% filter(dt_di_florio_cardiff$PHENOTYPE_FINE.x != dt_di_florio_cardiff$PHENOTYPE_FINE.y) %>% select(PARTICIPANT_ID, PHENOTYPE_FINE.y)
	where <- which(dt$PARTICIPANT_ID %in% samples_different$PARTICIPANT_ID)
	dt[where,]$PHENOTYPE_FINE <- samples_different$PHENOTYPE_FINE.y
}

# Send stupid values to NA
dt_di_florio_cardiff <- dt_di_florio_cardiff %>% 
	mutate(PSYCHOSIS.x = as.numeric(ifelse(PSYCHOSIS.x %in% c(0,1), PSYCHOSIS.x, NA)),
		   PSYCHOSIS.y = as.numeric(ifelse(PSYCHOSIS.y %in% c(0,1), PSYCHOSIS.y, NA)))

# Check agreement of psychosis where they are defined in both.
if (all(dt_di_florio_cardiff$PSYCHOSIS.x == dt_di_florio_cardiff$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

dt_di_florio_cardiff <- dt_di_florio_cardiff %>% mutate(PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x))

psychosis_cardiff <- merge(dt_di_florio_cardiff %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, PSYCHOSIS),
					 dt_cardiff %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, PSYCHOSIS), by = c("PARTICIPANT_ID", "PSYCHOSIS"), all=TRUE)
psychosis_cardiff <- psychosis_cardiff %>% mutate(LOCATION="Cardiff, UK")

# Check that this matches the psychosis data in the ICCBD folder.
dt_psych_ICCBD_header <- fread("ICCBD_data/bip_psychosis_v1.pheno.header", header=FALSE)
dt_psych_ICCBD <- fread("ICCBD_data/bip_psychosis_v1.pheno", header=FALSE)
names(dt_psych_ICCBD) <- unlist(dt_psych_ICCBD_header)
dt_psych_ICCBD <- dt_psych_ICCBD %>% rename(PARTICIPANT_ID = FID) %>% 
	mutate(PSYCHOSIS = ifelse(psychosis_sub==0, NA, ifelse(psychosis_sub==1, 0, 1))) %>% select(-psychosis_sub)

psychosis_cardiff <- merge(dt_psych_ICCBD, psychosis_cardiff, by='PARTICIPANT_ID', all=TRUE) %>% 
	mutate(study = ifelse(is.na(study), 'cardiff', study)) %>% 
	filter(study =='cardiff')

if (all(psychosis_cardiff$PSYCHOSIS.x == psychosis_cardiff$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

psychosis_cardiff <- psychosis_cardiff %>% mutate(PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x)) %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, LOCATION, PSYCHOSIS)
aao_cardiff <- dt %>% select(PARTICIPANT_ID, LOCATION, AGE_FI_24, AGE_FI_40, AGE_FS_24, AGE_FS_40, AGE_D_24, AGE_D_40)

# Swedish data.
cat("Clean and merge in Swedish information...\n\n")
dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv") %>% filter(LOCATION == "Stockholm, SWE")

# The following two files should contain the same data.
dt_ICCBD_stockholm <- fread("ICCBD_data/report_phenotype_sweden_25JUN2016.txt") %>% 
	select(
		PARTICIPANT_ID=participant_id, PSYCHOSIS=psychosis,
		# Adding in the age of onset variables as well
		# Age of first impairment, binned.
		AGE_FI_24 = age_FI_24, AGE_FI_40 = age_FI_40,
		# Age of first symptoms, binned.
		AGE_FS_24 = age_FS_24, AGE_FS_40 = Age_FS_40,
		# Age of first diagnosis, binned.
		AGE_D_24 = age_D_24, AGE_D_40 = age_D_40)

dt_smoller_stockholm_tmp <- fread("Smoller_data/lah1.output.pheno", colClasses="character")
dt_smoller_stockholm_tmp  <- dt_smoller_stockholm_tmp %>% 
	mutate(PARTICIPANT_ID=gsub(".*\\*", "", dt_smoller_stockholm_tmp$FID)) %>% 
	select(PARTICIPANT_ID, PSYCHOSIS=psychosis,
		# Adding in the age of onset variables as well
		# Age of first impairment, binned.
		AGE_FI_24 = age_fi_24, AGE_FI_40 = age_fi_40,
		# Age of first symptoms, binned.
		AGE_FS_24 = age_fs_24, AGE_FS_40 = age_fs_40,
		# Age of first diagnosis, binned.
		AGE_D_24 = age_d_24, AGE_D_40 = age_d_40)

swa2 <- fread("Smoller_data/swa2.output.pheno", colClasses="character") %>% 
		select(PARTICIPANT_ID=IID, PSYCHOSIS=psychosis,
			AGE_FI_24 = age_fi_24, AGE_FI_40 = age_fi_40,
			# Age of first symptoms, binned.
			AGE_FS_24 = age_fs_24, AGE_FS_40 = age_fs_40,
			# Age of first diagnosis, binned.
			AGE_D_24 = age_d_24, AGE_D_40 = age_d_40)

swei <- fread("Smoller_data/swei.output.pheno", colClasses="character") %>%
		select(PARTICIPANT_ID=IID, PSYCHOSIS=psychosis,
			AGE_FI_24 = age_fi_24, AGE_FI_40 = age_fi_40,
			# Age of first symptoms, binned.
			AGE_FS_24 = age_fs_24, AGE_FS_40 = age_fs_40,
			# Age of first diagnosis, binned.
			AGE_D_24 = age_d_24, AGE_D_40 = age_d_40)

dt_smoller_stockholm <- merge(dt_smoller_stockholm_tmp,
	merge(swa2, swei,
		by=c("PARTICIPANT_ID", "PSYCHOSIS",
			"AGE_FI_24", "AGE_FI_40", 
			"AGE_FS_24", "AGE_FS_40",
			"AGE_D_24", "AGE_D_40"), all=TRUE),
	by=c("PARTICIPANT_ID", "PSYCHOSIS",
		"AGE_FI_24", "AGE_FI_40",
		"AGE_FS_24", "AGE_FS_40",
		"AGE_D_24", "AGE_D_40"), all=TRUE)

dt_smoller_stockholm <- dt_smoller_stockholm %>% 
	mutate(PSYCHOSIS = as.numeric(PSYCHOSIS))
# Now, grab all the Stockholm information from these two merged files.
dt_stockholm <- merge(dt_ICCBD_stockholm, dt_smoller_stockholm, by="PARTICIPANT_ID", all=TRUE)
# Send stupid values to NA
dt_stockholm <- dt_stockholm %>% 
	mutate(PSYCHOSIS.x = ifelse(PSYCHOSIS.x %in% c(0,1), PSYCHOSIS.x, NA),
		   PSYCHOSIS.y = ifelse(PSYCHOSIS.y %in% c(0,1), PSYCHOSIS.y, NA),
		   AGE_FI_24.x = ifelse(AGE_FI_24.x %in% c(0,1,2), AGE_FI_24.x, NA),
		   AGE_FI_24.y = ifelse(AGE_FI_24.y %in% c(0,1,2), AGE_FI_24.y, NA),
		   AGE_FI_40.x = ifelse(AGE_FI_40.x %in% c(0,1,2), AGE_FI_40.x, NA),
		   AGE_FI_40.y = ifelse(AGE_FI_40.y %in% c(0,1,2), AGE_FI_40.y, NA),
		   AGE_FS_24.x = ifelse(AGE_FS_24.x %in% c(0,1,2), AGE_FS_24.x, NA),
		   AGE_FS_24.y = ifelse(AGE_FS_24.y %in% c(0,1,2), AGE_FS_24.y, NA),
		   AGE_FS_40.x = ifelse(AGE_FS_40.x %in% c(0,1,2), AGE_FS_40.x, NA),
		   AGE_FS_40.y = ifelse(AGE_FS_40.y %in% c(0,1,2), AGE_FS_40.y, NA),
		   AGE_D_24.x = ifelse(AGE_D_24.x %in% c(0,1,2), AGE_D_24.x, NA),
		   AGE_D_24.y = ifelse(AGE_D_24.y %in% c(0,1,2), AGE_D_24.y, NA),
		   AGE_D_40.x = ifelse(AGE_D_40.x %in% c(0,1,2), AGE_D_40.x, NA),
		   AGE_D_40.y = ifelse(AGE_D_40.y %in% c(0,1,2), AGE_D_40.y, NA))
		   
# Ensure that the psychosis information matches.
if (all(dt_stockholm$PSYCHOSIS.x == dt_stockholm$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

# Ensure that all the age of onset information matches.

# Impairment
if (all(dt_stockholm$AGE_FI_24.x == dt_stockholm$AGE_FI_24.y, na.rm=TRUE)) {
	cat("Age of onset; first impairment (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first impairment (24) does not agree in these two files\n")
}

if (all(dt_stockholm$AGE_FI_40.x == dt_stockholm$AGE_FI_40.y, na.rm=TRUE)) {
	cat("Age of onset; first impairment (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first impairment (40) does not agree in these two files\n")
}

# First symptoms
if (all(dt_stockholm$AGE_FS_24.x == dt_stockholm$AGE_FS_24.y, na.rm=TRUE)) {
	cat("Age of onset; first symptoms (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first symptoms (24) does not agree in these two files\n")
	print(dt_stockholm %>% filter(dt_stockholm$`AGE_FS_24.x` != dt_stockholm$`AGE_FS_24.y`))
	# Replace the discrepancy with the ICCBD version (this is the most up to date)
	where <- which(dt_stockholm$AGE_FS_24.x != dt_stockholm$AGE_FS_24.y)
	dt_stockholm$`AGE_FS_24.y`[where] <- dt_stockholm$`AGE_FS_24.x`[where]
}

# Check again.
if (all(dt_stockholm$AGE_FS_24.x == dt_stockholm$AGE_FS_24.y, na.rm=TRUE)) {
	cat("Age of onset; first symptoms (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first symptoms (24) does not agree in these two files\n")
}

if (all(dt_stockholm$AGE_FS_40.x == dt_stockholm$AGE_FS_40.y, na.rm=TRUE)) {
	cat("Age of onset; first symptoms (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; first symptoms (40) does not agree in these two files\n")
}

# Diagnosis
if (all(dt_stockholm$AGE_D_24.x == dt_stockholm$AGE_D_24.y, na.rm=TRUE)) {
	cat("Age of onset; diagnosis (24) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; diagnosis (24) does not agree in these two files\n")
}

if (all(dt_stockholm$AGE_D_40.x == dt_stockholm$AGE_D_40.y, na.rm=TRUE)) {
	cat("Age of onset; diagnosis (40) (where defined) agrees between these two files\n")
} else {
	cat("Age of onset; diagnosis (40) does not agree in these two files\n")
}

# Where one of the values is missing, pick the non-missing value
dt_stockholm <- dt_stockholm %>% 
	mutate(
		PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x),
		AGE_FI_24 = ifelse(is.na(AGE_FI_24.x), AGE_FI_24.y, AGE_FI_24.x), 
		AGE_FI_40 = ifelse(is.na(AGE_FI_40.x), AGE_FI_40.y, AGE_FI_40.x), 
		AGE_FS_24 = ifelse(is.na(AGE_FS_24.x), AGE_FS_24.y, AGE_FS_24.x), 
		AGE_FS_40 = ifelse(is.na(AGE_FS_40.x), AGE_FS_40.y, AGE_FS_40.x), 
		AGE_D_24 = ifelse(is.na(AGE_D_24.x), AGE_D_24.y, AGE_D_24.x), 
		AGE_D_40 = ifelse(is.na(AGE_D_40.x), AGE_D_40.y, AGE_D_40.x)) %>% 
	select(PARTICIPANT_ID, PSYCHOSIS,
		AGE_FI_24, AGE_FI_40, AGE_FS_24, AGE_FS_40, AGE_D_24, AGE_D_40)

# Check that the psychosis information that was sent over matches what we have (where there is an intersection).
dt_landen_psych <- fread("Landen_data/ICCBD_REPORT_PSYCHOSIS.tsv")

dt_stockholm <- merge(dt_stockholm, dt_landen_psych, by="PARTICIPANT_ID", all=TRUE) %>% 
	mutate(PSYCHOSIS.x = as.numeric(ifelse(PSYCHOSIS.x %in% c(0,1), PSYCHOSIS.x, NA)),
		   PSYCHOSIS.y = as.numeric(ifelse(PSYCHOSIS.y %in% c(0,1), PSYCHOSIS.y, NA)))

if (all(dt_stockholm$PSYCHOSIS.x == dt_stockholm$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

dt_stockholm <- dt_stockholm %>% mutate(PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x))

psychosis_stockholm <- dt_stockholm %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, PSYCHOSIS)
psychosis_stockholm <- psychosis_stockholm %>% mutate(LOCATION="Stockholm, SWE")
aao_stockholm <- dt_stockholm %>% select(PARTICIPANT_ID, AGE_FI_24, AGE_FI_40, AGE_FS_24, AGE_FS_40, AGE_D_24, AGE_D_40)
aao_stockholm <- aao_stockholm %>% mutate(LOCATION="Stockholm, SWE")
# Check that this matches the psychosis data in the ICCBD folder. 
dt_psych_ICCBD_header <- fread("ICCBD_data/bip_psychosis_v1.pheno.header", header=FALSE)
dt_psych_ICCBD <- fread("ICCBD_data/bip_psychosis_v1.pheno", header=FALSE)
names(dt_psych_ICCBD) <- unlist(dt_psych_ICCBD_header)
dt_psych_ICCBD <- dt_psych_ICCBD %>% rename(PARTICIPANT_ID = FID) %>% 
	mutate(PSYCHOSIS = ifelse(psychosis_sub==0, NA, ifelse(psychosis_sub==1, 0, 1))) %>% select(-psychosis_sub)

psychosis_stockholm <- merge(dt_psych_ICCBD, psychosis_stockholm, by='PARTICIPANT_ID', all=TRUE) %>% 
	mutate(study = ifelse(is.na(study), 'ki', study)) %>% 
	filter(study =='ki')

if (all(dt_stockholm$PSYCHOSIS.x == dt_stockholm$PSYCHOSIS.y, na.rm=TRUE)) {
	cat("Psychosis information (where defined) agrees between these two files\n")
} else {
	cat("Psychosis information does not agree in these two files\n")
}

psychosis_stockholm <- psychosis_stockholm %>% mutate(PSYCHOSIS = ifelse(is.na(PSYCHOSIS.x), PSYCHOSIS.y, PSYCHOSIS.x)) %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, LOCATION, PSYCHOSIS)

# MGH phenotype data.
cat("Clean and merge in Boston information...\n\n")
dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv") %>% filter(LOCATION == "Boston, USA")

dt_smoller_mgh <- fread("ICCBD_data/ICCBD.Phenotype.MGH.02Apr14_edit.txt") %>%
	mutate(ICCBD_PARTICIPANT_ID = as.character(ParticipantID), PSYCHOSIS=Psychosis,
		PHENOTYPE_COARSE = case_when(cohort==0 ~ "Control", 
			SCZ_BP == 1 ~ "Schizoaffective",
			cohort == 1 ~ "Bipolar Disorder"),
		PHENOTYPE_FINE = case_when(
			cohort == 0 ~ "Control",
			BP1 == 1 ~ "Bipolar Disorder 1",
			BP2 == 1 ~ "Bipolar Disorder 2",
			BP_NOS == 1 ~ "Bipolar Disorder NOS",
			SCZ_BP == 1 ~ "Schizoaffective",
			cohort == 1 ~ "Bipolar Disorder")
		)

dt_key <- fread("ICCBD_data/iccbd_full_manifest_final.csv") %>% rename(SAMPLE_ALIAS = CrimsonSampleID, ICCBD_PARTICIPANT_ID=patient_num) %>% select(SAMPLE_ALIAS, ICCBD_PARTICIPANT_ID)

# Need to merge in the PT-ID or 
dt_smoller_mgh <- merge(dt_smoller_mgh, dt_key, by="ICCBD_PARTICIPANT_ID")

dim(dt)
dt_boston <- merge((dt %>% mutate(SAMPLE_ALIAS = gsub("_A", "", SAMPLE_ALIAS))), dt_smoller_mgh, by="SAMPLE_ALIAS", all.x=TRUE)

# Where things are not missing, ensure that the phenotypes match.
if (all(dt_boston$PHENOTYPE_COARSE.x == dt_boston$PHENOTYPE_COARSE.y, na.rm=TRUE)) {
	cat("Coarse phenotype information (where defined) agrees between these two files\n")
} else {
	cat("Coarse phenotype information does not agree in these two files\n")
}

if (all(dt_boston$PHENOTYPE_FINE.x == dt_boston$PHENOTYPE_FINE.y, na.rm=TRUE)) {
	cat("Fine phenotype information (where defined) agrees between these two files\n")
} else {
	cat("Fine phenotype information does not agree in these two files\n")
	print(dt_boston %>% filter(PHENOTYPE_FINE.x != PHENOTYPE_FINE.y) %>% select(PARTICIPANT_ID, PHENOTYPE_FINE.x, PHENOTYPE_FINE.y))
	# The new data has some extra information!
	dt_boston[which(dt_boston$PHENOTYPE_FINE.x != dt_boston$PHENOTYPE_FINE.y),]$PHENOTYPE_FINE.x <- dt_boston[which(dt_boston$PHENOTYPE_FINE.x != dt_boston$PHENOTYPE_FINE.y),]$PHENOTYPE_FINE.y
	# So we just need to read in the file and change this single entry.
	dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv")
	dt[dt$PARTICIPANT_ID == dt_boston[which(dt_boston$PHENOTYPE_FINE.x != dt_boston$PHENOTYPE_FINE.y),]$PARTICIPANT_ID,]$PHENOTYPE_FINE <- dt_boston[which(dt_boston$PHENOTYPE_FINE.x != dt_boston$PHENOTYPE_FINE.y),]$PHENOTYPE_FINE.y
	fwrite(dt, file="BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv", row.names=FALSE, sep='\t')
}

# Check where things are missing.
cat("Checking where things are missing...\n")
dt_boston_missing <- dt %>% filter(SAMPLE_ALIAS %in% setdiff(gsub("_A", "", dt$SAMPLE_ALIAS), dt_smoller_mgh$SAMPLE_ALIAS))
dt_boston_missing <- merge(dt_boston_missing, dt_key, by="SAMPLE_ALIAS") %>% select(ICCBD_PARTICIPANT_ID)

# Some are in the usa2 pheno file, others are in the i2b2 fam file, but not present in the pheno file.

all_fams <- merge(
	merge(
		merge(
			fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.fam", colClasses="character") %>% mutate(file="i2b2"),
			fread("Smoller_data/found_on_cluster/smol4v11_pchip_updateid.gencall.fam", colClasses="character") %>% mutate(file="smol4v11"),
			all=TRUE),
		merge(
			fread("Smoller_data/found_on_cluster/usa1s_pchip_updateid.gencall.fam", colClasses="character") %>% mutate(file="usa1"),
			fread("Smoller_data/found_on_cluster/usa2smoll_pchip_updateid.gencall.fam", colClasses="character") %>% mutate(file="usa2"),
			all=TRUE), all=TRUE
		),
	fread("Smoller_data/found_on_cluster/usa3con_pchip_updateid.gencall.fam", colClasses="character") %>% mutate(file="usa3"), all=TRUE
)

all_fams %>% filter(all_fams$V2 %in% dt_boston_missing$ICCBD_PARTICIPANT_ID)

# Are they in the phenotype files?

dt_pheno <- merge(fread("Smoller_data/found_on_cluster/usa2smoll_pchip_updateid.gencall.pheno") %>% mutate(file="usa2"), fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.pheno") %>% mutate(file="i2b2"), all=TRUE)
names(dt_pheno)<- c("PTID", "ParticipantID", "famgender", "race", "bip_NLP_single", "bip_Strict_single", "bip_Broad_single", "bip_SV_single", "bip_NLP_nohier", "bip_Strict_nohier", "bip_Broad_nohier", "bip_SV_nohier", "bip_NLP_hier", "bip_Strict_hier", "bip_Broad_hier", "bip_SV_hier", "BP1", "BP2", "BPNOS", "BPSCZ", "bip_hier_final", "BPTYPE_final", "Panic", "Anxiety", "Sub_dep", "Alc_dep", "ADHD", "ASD", "Dev_Del", "Int_Dis", "Epilepsy", "Migraines", "Psychosis", "Suicidality", "Age_FI_24", "Age_FI_40", "Age_FS_24", "Age_FS_40", "Age_D_24", "Age_D_40", "FHX_BP", "FHX_Mood", "FHX_Psy", "LI", "Delivery_num", "PostP_Psychosis", "PostP_Dep", "Peri_mood", "Neurocog", "file")

dt_boston_missing$ICCBD_PARTICIPANT_ID %in% dt_pheno$ParticipantID
dt_pheno %>% filter(dt_pheno$ParticipantID %in% dt_boston_missing$ICCBD_PARTICIPANT_ID)
# These also have the phenotype information missing. Good. Nothing left to add.

setdiff((all_fams %>% filter(file=="i2b2"))$V2, fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.pheno")$V2) %in% dt_boston_missing$ICCBD_PARTICIPANT_ID
# Are they 'missing phenotype information' - Yes!
fread("Smoller_data/found_on_cluster/i2b2_pchip_updateid.gencall.missing.pheno.txt")$ParticipantID %in% dt_boston_missing$ICCBD_PARTICIPANT_ID

# So leave these as is. Though a large proportion of them are classified as cases on BSP, we don't have phenotype information for them.
# Psychosis information that we wish to merge in is:
cat('Obtain the psychosis information from the Boston cohorts...\n')
psychosis_boston <- dt_boston %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, LOCATION, PSYCHOSIS)
aao_boston <- dt_boston %>% select(
	PARTICIPANT_ID,
	LOCATION,
	AGE_FI_24 = Age_FI_24,
	AGE_FI_40 = Age_FI_40,
	AGE_FS_24 = Age_FS_24,
	AGE_FS_40 = Age_FS_40,
	AGE_D_24 = Age_D_24,
	AGE_D_40 = Age_D_40
)
# Finally, remove the duplicates. Exactly one of these actually matters, all the rest are controls or don't 
# have AAO information - for each, pick the first one.
duplicated_PARTICIPANT_IDs <- aao_boston$PARTICIPANT_ID[which(duplicated(aao_boston$PARTICIPANT_ID))]
for(id in duplicated_PARTICIPANT_IDs) {
	where <- which(aao_boston$PARTICIPANT_ID == id)
	aao_boston <- aao_boston[-where[2:length(where)],]
}

# German data.
cat("Clean and merge in Cardiff information...\n\n")
psychosis_german <- fread("Reif_data/Reif_BD_Info_toDuncan_672019_TMK.txt") %>% rename(PSYCHOSIS = `Psychosis (Y/N)`) %>% select(SAMPLE_ALIAS, PSYCHOSIS) %>% filter(!is.na(PSYCHOSIS))
psychosis_german <- merge(psychosis_german %>% mutate(PSYCHOSIS = ifelse(PSYCHOSIS == "Y", 1, 0)), dt, all.x=TRUE) %>% select(PARTICIPANT_ID, LOCATION, PSYCHOSIS)

# Andrew McQuillin data. Using the OPCRIT information - This is incorrect, need to use the older version of the data for this.

# dt_mcquillin_london <- fread("McQuillin_data/annabel_121119_04pm.pheno") %>% rename(SAMPLE_ALIAS=id)
# dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv")
# dt_london <- merge(dt_mcquillin_london, dt %>% filter(PI=="Andrew McQuillin"), by="SAMPLE_ALIAS", all.y=TRUE)

# # Check
# dt_london %>% select(PHENOTYPE_COARSE, PHENOTYPE_FINE, bp_type_inferred) %>%
# 	filter(!is.na(bp_type_inferred)) %>% 
# 	mutate(PHENOTYPE_FINE_CHECK = case_when(
# 		bp_type_inferred == "BP1" ~ "Bipolar Disorder 1",
# 		bp_type_inferred == "bp1" ~ "Bipolar Disorder 1",
# 		bp_type_inferred == "BP2" ~ "Bipolar Disorder 2",
# 		bp_type_inferred == "bp2" ~ "Bipolar Disorder 2",
# 		bp_type_inferred == "BP9" ~ "Bipolar Disorder",
# 		bp_type_inferred == "SABP" ~ "Schizoaffective",
# 		bp_type_inferred == "SADEP" ~ "Schizoaffective")
# 	)


# # Check agreement of phenotypes and subtypes where they are defined.
# if (all(dt_london$PHENOTYPE_FINE == dt_london$PHENOTYPE_FINE_CHECK, na.rm=TRUE)) {
# 	cat("Fine phenotype information (where defined) agrees between these two files\n")
# } else {
# 	cat("Fine phenotype information does not agree in these two files\n")	
# }

# dt_london <- dt_london %>% mutate(PSYCHOSIS =
# 	ifelse(
# 		(is.na(opcrit.52) & is.na(opcrit.54) & is.na(opcrit.55) & is.na(opcrit.57) & is.na(opcrit.58) &
# 		is.na(opcrit.59) & is.na(opcrit.60) & is.na(opcrit.61) & is.na(opcrit.62) & is.na(opcrit.63) &
# 		is.na(opcrit.64) & is.na(opcrit.65) & is.na(opcrit.66) & is.na(opcrit.67) & is.na(opcrit.68) &
# 		is.na(opcrit.69) & is.na(opcrit.70) & is.na(opcrit.71) & is.na(opcrit.72) & is.na(opcrit.73) &
# 		is.na(opcrit.74) & is.na(opcrit.75) & is.na(opcrit.76) & is.na(opcrit.77)
# 		), NA, 0)
# 	)

# psychosis_london <- dt_london %>% mutate(PSYCHOSIS_SCALE = PSYCHOSIS + ((opcrit.52 %in% c(1,2,3)) + (opcrit.54 == 1) + 
# 	(opcrit.55 == 1) + (opcrit.57 %in% c(1,2)) + (opcrit.58 == 1) + (opcrit.59 == 1) + (opcrit.60 == 1) + (opcrit.61 == 1) +
# 	(opcrit.62 == 1) + (opcrit.63 == 1) + (opcrit.64 == 1) + (opcrit.65 == 1) + (opcrit.66 == 1) + (opcrit.67 == 1) +
# 	(opcrit.68 == 1) + (opcrit.69 == 1) + (opcrit.70 == 1) + (opcrit.71 == 1) + (opcrit.72 == 1) + (opcrit.73 == 1) + 
# 	(opcrit.74 == 1) + (opcrit.75 == 1) + (opcrit.76 == 1) + (opcrit.77 == 1))) %>% mutate(PSYCHOSIS = ifelse(PSYCHOSIS_SCALE > 0, 1, 0)) %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, PSYCHOSIS)

# Fixing the McQuillin data - check 01_checking_mcqullin_data - here I harmonise the two older datasets.
cat("Clean and merge in London information...\n\n")
source("01_checking_mcquillin_data.r")
dt_mcquillin_london <- fread("McQuillin_data/harmonised_opcrit_mcquillin_july_sept_nov.tsv") %>% rename(SAMPLE_ALIAS=id)
dt <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv")
dt_london <- merge(dt_mcquillin_london, dt %>% filter(PI=="Andrew McQuillin"), by="SAMPLE_ALIAS", all.y=TRUE)

# Check
dt_london %>% select(PHENOTYPE_COARSE, PHENOTYPE_FINE, bp_type_inferred) %>%
	filter(!is.na(bp_type_inferred)) %>% 
	mutate(PHENOTYPE_FINE_CHECK = case_when(
		bp_type_inferred == "BP1" ~ "Bipolar Disorder 1",
		bp_type_inferred == "bp1" ~ "Bipolar Disorder 1",
		bp_type_inferred == "BP2" ~ "Bipolar Disorder 2",
		bp_type_inferred == "bp2" ~ "Bipolar Disorder 2",
		bp_type_inferred == "BP9" ~ "Bipolar Disorder",
		bp_type_inferred == "SABP" ~ "Schizoaffective",
		bp_type_inferred == "SADEP" ~ "Schizoaffective")
)

# Check agreement of phenotypes and subtypes where they are defined.
if (all(dt_london$PHENOTYPE_FINE == dt_london$PHENOTYPE_FINE_CHECK, na.rm=TRUE)) {
	cat("Fine phenotype information (where defined) agrees between these two files\n")
} else {
	cat("Fine phenotype information does not agree in these two files\n")	
}

dt_london <- dt_london %>% mutate(PSYCHOSIS =
	ifelse(
		(is.na(opcrit.52) & is.na(opcrit.54) & is.na(opcrit.55) & is.na(opcrit.57) & is.na(opcrit.58) &
		is.na(opcrit.59) & is.na(opcrit.60) & is.na(opcrit.61) & is.na(opcrit.62) & is.na(opcrit.63) &
		is.na(opcrit.64) & is.na(opcrit.65) & is.na(opcrit.66) & is.na(opcrit.67) & is.na(opcrit.68) &
		is.na(opcrit.69) & is.na(opcrit.70) & is.na(opcrit.71) & is.na(opcrit.72) & is.na(opcrit.73) &
		is.na(opcrit.74) & is.na(opcrit.75) & is.na(opcrit.76) & is.na(opcrit.77)
		), NA, 0)
	)

psychosis_london <- dt_london %>% mutate(PSYCHOSIS_SCALE = PSYCHOSIS + ((opcrit.52 %in% c(1,2,3)) + (opcrit.54 == 1) + 
	(opcrit.55 == 1) + (opcrit.57 %in% c(1,2)) + (opcrit.58 == 1) + (opcrit.59 == 1) + (opcrit.60 == 1) + (opcrit.61 == 1) +
	(opcrit.62 == 1) + (opcrit.63 == 1) + (opcrit.64 == 1) + (opcrit.65 == 1) + (opcrit.66 == 1) + (opcrit.67 == 1) +
	(opcrit.68 == 1) + (opcrit.69 == 1) + (opcrit.70 == 1) + (opcrit.71 == 1) + (opcrit.72 == 1) + (opcrit.73 == 1) + 
	(opcrit.74 == 1) + (opcrit.75 == 1) + (opcrit.76 == 1) + (opcrit.77 == 1))) %>% mutate(PSYCHOSIS = ifelse(PSYCHOSIS_SCALE > 0, 1, 0)) %>% filter(!is.na(PSYCHOSIS)) %>% select(PARTICIPANT_ID, LOCATION, PSYCHOSIS)

aao_london <- dt_london %>% mutate(
	AGE_FI_24 = 
		case_when(
			opcrit.04 < 12 ~ 0,
			((opcrit.04 >= 12) & (opcrit.04 <= 24)) ~ 1,
			(opcrit.04 > 24) ~ 2
		),
	AGE_FI_40 = 
		case_when(
			opcrit.04 < 18 ~ 0,
			((opcrit.04 >= 18) & (opcrit.04 <= 40)) ~ 1,
			(opcrit.04 > 40) ~ 2
		),
	AGE_FS_24 = NA,
	AGE_FS_40 = NA,
	AGE_D_24 = NA,
	AGE_D_40 = NA) %>% select(PARTICIPANT_ID, LOCATION, AGE_FI_24, AGE_FI_40, AGE_FS_24, AGE_FS_40, AGE_D_24, AGE_D_40)

psychosis <- merge(merge(merge(psychosis_cardiff, psychosis_stockholm, all=TRUE), merge(psychosis_boston, psychosis_german, all=TRUE), all=TRUE), psychosis_london, all=TRUE)
aao <- merge(merge(aao_cardiff, aao_stockholm, all=TRUE), merge(aao_london, aao_boston, all=TRUE), all=TRUE)
# # Last piece - I've found a "combined phenotype file on the cluster" - let's check that this agrees with what we have in this matrix (where it matches).
dt <- merge(dt, psychosis, by=c("PARTICIPANT_ID","LOCATION"), all.x=TRUE)
dt <- merge(dt, aao, by=c("PARTICIPANT_ID","LOCATION"), all.x=TRUE)

dt_key <- fread("ICCBD_data/iccbd_full_manifest_final.csv") %>% rename(SAMPLE_ALIAS = CrimsonSampleID, ICCBD_PARTICIPANT_ID=patient_num) %>% select(SAMPLE_ALIAS, ICCBD_PARTICIPANT_ID)

combined_ICCBD <- fread("ICCBD_data/phenotype_merge_v4.txt") %>% mutate(PHENOTYPE_COARSE = ifelse(type == 0, "Control", "Bipolar Disorder"),
	PHENOTYPE_FINE = case_when(
		BP1==1 ~ "Bipolar Disorder 1",
		BP2==1 ~ "Bipolar Disorder 2",
		BP_NOS==1 ~ "Bipolar Disorder NOS",
		SCZ_BP==1 ~ "Schizoaffective",
		type==0 ~ "Control",
		type==1 ~ "Bipolar Disorder")) %>% rename(PARTICIPANT_ID = participant_id)
combined_ICCBD$PHENOTYPE_COARSE[combined_ICCBD$PHENOTYPE_FINE == "Schizoaffective"] <- "Schizoaffective"

# Need to merge in the PT-ID or 
dt_check <- merge(dt, dt_key, by="SAMPLE_ALIAS", all.x=TRUE)
dt_check <- dt_check %>% mutate(new_key = ifelse(ICCBD_PARTICIPANT_ID %in% combined_ICCBD$PARTICIPANT_ID, ICCBD_PARTICIPANT_ID, PARTICIPANT_ID))

# The differences are the Karolinska samples that underwent rediagnosis.
dt_ICCBD_check <- merge(combined_ICCBD %>% rename(new_key = PARTICIPANT_ID), dt_check, by="new_key", all.y=TRUE)
dt_ICCBD_check %>% select(PHENOTYPE_COARSE.x, PHENOTYPE_COARSE.y, study) %>% filter(PHENOTYPE_COARSE.x != PHENOTYPE_COARSE.y)
dt_ICCBD_check %>% select(PHENOTYPE_FINE.x, PHENOTYPE_FINE.y, study) %>% filter(PHENOTYPE_FINE.x != PHENOTYPE_FINE.y)
dt_ICCBD_check %>% select(PARTICIPANT_ID, PSYCHOSIS, psychosis, study) %>% filter(PSYCHOSIS != psychosis)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_FI_24, age_FI_24, study) %>% filter(AGE_FI_24 != age_FI_24)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_FI_40, age_FI_40, study) %>% filter(AGE_FI_40 != age_FI_40)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_FS_24, age_FS_24, study) %>% filter(AGE_FS_24 != age_FS_24)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_FS_40, age_FS_40, study) %>% filter(AGE_FS_40 != age_FS_40)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_D_24, age_D_24, study) %>% filter(AGE_D_24 != age_D_24)
dt_ICCBD_check %>% select(PARTICIPANT_ID, AGE_D_40, age_D_40, study) %>% filter(AGE_D_40 != age_D_40)

# Psychosis check
all(is.na(dt_ICCBD_check$PSYCHOSIS) == is.na(dt_ICCBD_check$psychosis))
which(dt_ICCBD_check$PSYCHOSIS != dt_ICCBD_check$psychosis)

# AAO check
aao_check <- function(dt_ICCBD_check, variable_dt, variable_ICCBD_dt) {

	cat(all(is.na(dt_ICCBD_check[[variable_dt]]) == is.na(dt_ICCBD_check[[variable_ICCBD_dt]])),'\n')
	# print(dt_ICCBD_check %>% filter(is.na(dt_ICCBD_check[[variable_dt]]) != is.na(dt_ICCBD_check[[variable_ICCBD_dt]])) %>% select(PARTICIPANT_ID, LOCATION, !!variable_dt, !!variable_ICCBD_dt))
	where_diff <- which(dt_ICCBD_check[[variable_dt]] != dt_ICCBD_check[[variable_ICCBD_dt]])
	cat('Differences at non-missing sample data:', length(where_diff), '\n')
	print(dt_ICCBD_check[where_diff,] %>% select(PARTICIPANT_ID, LOCATION, PI, !!variable_dt, !!variable_ICCBD_dt))
	
	where <- (is.na(dt_ICCBD_check[[variable_dt]]) & (!is.na(dt_ICCBD_check[[variable_ICCBD_dt]])))
	if(any(where)) {
		cat('Number of AAO gained:', sum(where), '\n')
		dt_ICCBD_check[[variable_dt]][where] <- dt_ICCBD_check[[variable_ICCBD_dt]][where]
	} else {
		cat('No AAO gained!\n')
	}
	return(dt_ICCBD_check)
}

# Let's merge in the extra information from the ICCBD combined file and report the extra information as we go.

# First impairment
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_FI_24', 'age_FI_24')
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_FI_40', 'age_FI_40')
# First symptoms
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_FS_24', 'age_FS_24')
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_FS_40', 'age_FS_40')
# Diagnosis
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_D_24', 'age_D_24')
dt_ICCBD_check <- aao_check(dt_ICCBD_check, 'AGE_D_40', 'age_D_40')

# Nothing extra to add for subphenotype information
options(width=150)
check <- aao_check(dt_ICCBD_check, 'PHENOTYPE_FINE.y', 'PHENOTYPE_FINE.x')
check <- aao_check(dt_ICCBD_check, 'PHENOTYPE_COARSE.y', 'PHENOTYPE_COARSE.x')

# Good - they all match where data is present (except for at tiny subset that we already know about), and we
# have gained a little more information from the phenotype_merge_v4 file.

dt_x <- merge(combined_ICCBD %>% rename(new_key = PARTICIPANT_ID), dt_check, by="new_key", all.x=TRUE)
dt_y <- merge(combined_ICCBD %>% rename(new_key = PARTICIPANT_ID), dt_check, by="new_key", all.y=TRUE)

dt_x %>% filter(new_key %in% setdiff(dt_x$new_key, dt_y$new_key)) %>% select(study, PHENOTYPE_FINE.x, PHENOTYPE_COARSE.x)
dt_y %>% filter(new_key %in% setdiff(dt_y$new_key, dt_x$new_key)) %>% select(LOCATION, PI, PHENOTYPE_FINE.y, PHENOTYPE_COARSE.y)

# Everything checked. May need to rethink the Swedish portion at a later date.

# Double check that everything that is Schizoaffective in fine is Shizoaffective in coarse and vice versa
table(dt %>% filter(PHENOTYPE_FINE == "Schizoaffective") %>% select(PHENOTYPE_COARSE))
table(dt %>% filter(PHENOTYPE_COARSE == "Schizoaffective") %>% select(PHENOTYPE_FINE))

dt <- dt_ICCBD_check %>% 
	mutate(
		PHENOTYPE_COARSE = PHENOTYPE_COARSE.y,
		PHENOTYPE_FINE = PHENOTYPE_FINE.y,
		PSYCHOSIS = as.logical(as.integer(PSYCHOSIS)),
		AGE_FI_24 = as.integer(AGE_FI_24),
		AGE_FI_40 = as.integer(AGE_FI_40),
		AGE_FS_24 = as.integer(AGE_FS_24),
		AGE_FS_40 = as.integer(AGE_FS_40),
		AGE_D_24 = as.integer(AGE_D_24),
		AGE_D_40 = as.integer(AGE_D_40)
	) %>%
	select(
		PARTICIPANT_ID,
		SAMPLE_ALIAS,
		GENDER,
		PROJECT_OR_COHORT,
		VCF,
		BAM,
		PCT_CONTAMINATION,
		PCT_CHIMERAS,
		LOCATION,
		PI, 
		INSTITUTION,
		PHENOTYPE_COARSE,
		PHENOTYPE_FINE,
		PSYCHOSIS,
		AGE_FI_24,
		AGE_FI_40,
		AGE_FS_24,
		AGE_FS_40,
		AGE_D_24,
		AGE_D_40
	)

# Apart from that, write the merged dt_check to disk.
fwrite(dt, file="BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv", row.names=FALSE, sep='\t', na='NA')

# Double check against old version
dt_1 <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv")
dt_2 <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_final.tsv")

print(dim(dt_1))
print(dim(dt_2))


