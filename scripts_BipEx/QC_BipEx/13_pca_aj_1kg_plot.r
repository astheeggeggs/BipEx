library(ggplot2)
library(dplyr)
library(ggsci)
library(data.table)
library(randomForest)

source('r_functions_and_parameters/pretty_plotting.r')
source("r_functions_and_parameters/r_options_BipEx.r")

# Throughout, we only consider the strictly defined European set.
save_figures <- TRUE

URV_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv'

df_URVs <- fread(URV_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df_URVs)[names(df_URVs) == 'Samples'] <- 's'

# The European samples with the Ashkenasi Jews, and 1kg. 
df_1kg_EUR_AJ <- fread(PCA_EUR_1KG_AJ_SCORES, sep='\t', header=TRUE, data.table=FALSE)
df_1kg_EUR_AJ <- df_1kg_EUR_AJ[sample(nrow(df_1kg_EUR_AJ), replace=FALSE),]
df_1kg_EUR_AJ <- merge(df_1kg_EUR_AJ, df_URVs, all.x=TRUE)

df_AJ <- fread('../../samples_BipEx/12_aj.sample_list', sep='\t', data.table=FALSE, header=FALSE)
df_1kg_EUR_AJ$POPULATION[df_1kg_EUR_AJ$s %in% df_AJ$V1] <- 'AJ'
df_1kg_EUR_AJ <- mutate(df_1kg_EUR_AJ, POPULATION = factor(POPULATION))

# All the European samples. 
df_1kg_EUR <- fread(PCA_EUR_1KG_SCORES, sep='\t', header=TRUE, data.table=FALSE)
df_1kg_EUR <- df_1kg_EUR[sample(nrow(df_1kg_EUR), replace=FALSE),]
df_1kg_EUR <- merge(df_1kg_EUR, df_URVs, all.x=TRUE)

# All strictly defined Europeans.

# Dataframes
dataframes <- c("df_1kg_EUR_AJ",
                "df_1kg_EUR")
filenames <- c("EUR_and_AJ_1kg",
               "EUR_1kg")

PCs <- c(1,3,5)

for(j in 1:length(dataframes)) {
    for (i in PCs) {
        # Colour by 1kg population, and case status.
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='POPULATION')
        p <- create_pretty_scatter(eval(parse(text = dataframes[j])), aes,
          add_final_layer=TRUE,
          final_layer=subset(eval(parse(text = dataframes[j])), SUPER_POPULATION == 'EUR' | POPULATION == 'AJ'),
          save_figure=save_figures, file=paste0(PLOTS, "13_PC", i, "PC", i+1, "_", filenames[j]),
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1))

        # Colour by sequencing collection.
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='LOCATION')
        p <- create_pretty_scatter(eval(parse(text = dataframes[j])), aes,
          add_final_layer=TRUE,
          final_layer=subset(eval(parse(text = dataframes[j])), SUPER_POPULATION == 'EUR'),
          save_figure=save_figures, file=paste0(PLOTS, "13_PC", i, "PC", i+1, "_", filenames[j], "_collection"),
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1))

        # Colour by sequencing collection.
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='PHENOTYPE_COARSE')
        p <- create_pretty_scatter(eval(parse(text = dataframes[j])), aes,
          add_final_layer=TRUE,
          final_layer=subset(eval(parse(text = dataframes[j])), SUPER_POPULATION == 'EUR'),
          save_figure=save_figures, file=paste0(PLOTS, "13_PC", i, "PC", i+1, "_", filenames[j], "_case_control"),
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1))

    }
}

df_train <- filter(df_1kg_EUR_AJ, (PHENOTYPE_COARSE == '1KG' | POPULATION == 'AJ')) %>%
    mutate(AJ = factor(POPULATION == 'AJ', labels=c('non-AJ', 'AJ'))) %>%
    select(c(PHENOTYPE_COARSE, AJ, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

df_predict <- filter(df_1kg_EUR_AJ, PHENOTYPE_COARSE != '1KG' & POPULATION != 'AJ') %>%
    select(c(s, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,6))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=df_train[PCs_to_use], y=df_train$AJ, ntree=5000)
rf_probs <- predict(rf, df_predict[PCs_to_use], type='prob')

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, 0.95)
classification <- as.character(predict(rf, df_predict[PCs_to_use]))
df_predict$classification_loose <- factor(classification)
classification[unsure] <- 'unsure'
df_predict$classification_strict <- factor(classification)

df_classify_AJ <- df_predict %>% 
    select(s, classification_strict, classification_loose) %>% 
    full_join(df_1kg_EUR_AJ, by='s') %>% filter(PHENOTYPE_COARSE != '1KG')

# The number of non-European samples
print(paste0("Number of AJ: ", table((df_classify_AJ %>% filter(classification_strict == 'AJ'))$PHENOTYPE_FINE)))
print("Number of Unsure:\n")
print(table((df_classify_AJ %>% filter(classification_strict == 'unsure'))$PHENOTYPE_FINE))

df_classify_AJ$classification_loose <- as.character(df_classify_AJ$classification_loose)
df_classify_AJ$classification_strict <- as.character(df_classify_AJ$classification_strict)
df_classify_AJ$classification_loose[which(is.na(df_classify_AJ$classification_loose))] <- "AJ trained"
df_classify_AJ$classification_strict[which(is.na(df_classify_AJ$classification_strict))] <- "AJ trained"
df_classify_AJ$classification_loose <- as.factor(df_classify_AJ$classification_loose)
df_classify_AJ$classification_strict <- as.factor(df_classify_AJ$classification_strict)

for (i in PCs) {
    aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
    create_pretty_scatter(df_classify_AJ, aes, save_figure=save_figures,
        file=paste0(PLOTS,'13_PC',i,'_PC',i+1,'_classify_AJ'),
        x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
}

fwrite(data.frame((df_classify_AJ %>% filter(classification_strict == 'AJ'))$s),
  file = "../../samples_BipEx/13_aj_classify.sample_list", col.names=FALSE)

