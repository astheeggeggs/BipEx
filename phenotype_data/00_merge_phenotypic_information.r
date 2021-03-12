library(data.table)
library(dplyr)

# Merge together all of the BSP files
df <- fread('BSP_phenotype_data_1.txt', na.strings = c("", "NA"), colClasses="character")
df2 <- fread('dalio_39k_from_bsp.txt', na.strings = c("", "NA"), colClasses="character")
where_dups <- which(duplicated(names(df2)))
names(df2)[where_dups] <- paste0(names(df2)[where_dups], '_2')
df <- merge(df, df2, all=TRUE)
df <- unique(df)

names(df)[names(df) == "Collection"] <- "COLLECTION"
names(df)[names(df) == "Primary Disease"] <- "PRIMARY_DISEASE"
names(df)[names(df) == "Gender"] <- "GENDER"
names(df)[which(names(df) == 'Sample ID')] <- 'SAMPLE_ID'
names(df)[which(names(df) == 'Collaborator Sample ID')] <- 'SAMPLE_ALIAS'
names(df)[which(names(df) == 'Participant ID(s)')] <- 'PARTICIPANT_ID'

df <- df %>% select("COLLECTION", "SAMPLE_ALIAS", "PRIMARY_DISEASE", "GENDER", "SAMPLE_ID", "PARTICIPANT_ID")
setkey(df, 'SAMPLE_ID')

df_key <- fread('DSP_data.txt', na.strings = c("", "NA"))
# Remove all rows that have an empty field for the SAMPLE_ID (SMID)
# These were repeated, because not sufficient aliquot.
df_key <- df_key %>% filter(!is.na(SAMPLE_ID))
setkey(data.table(df_key), 'SAMPLE_ID')

# Merge the result together
df <- merge(df, df_key, all=TRUE)

# Now remove duplicate rows after removal of SAMPLE_ID.
df <- unique(df %>% select(-SAMPLE_ID))

fwrite(df, file="BIP_phenotype_information.tsv", sep='\t')

# Got these files through the alternative means of contacting BSP directly.

df_check1 <- fread('PO-17611BspPart1.txt', na.strings = c("", "NA"))
df_check2 <- fread('PO-17611BspPart2.txt', na.strings = c("", "NA"))

df_check <- merge(df_check1, df_check2, all=TRUE)
df_check <- unique(df_check)

names(df_check)[names(df_check) == "Collection"] <- "COLLECTION"
names(df_check)[names(df_check) == "Primary Disease"] <- "PRIMARY_DISEASE"
names(df_check)[names(df_check) == "Gender"] <- "GENDER"
names(df_check)[which(names(df_check) == 'Sample ID')] <- 'SAMPLE_ID'
names(df_check)[which(names(df_check) == 'Collaborator Sample ID')] <- 'SAMPLE_ALIAS'

# The resultant file is the same! Great!