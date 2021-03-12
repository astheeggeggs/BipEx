library(ggplot2)
library(dplyr)
library(ggsci)
library(gridExtra)
library(randomForest)
library(data.table)

# Load plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')

# Get file locations, plotting locations, and thresholds.
source("r_functions_and_parameters/r_options_BipEx.r")

save_figures <- TRUE
perform_plotting <- TRUE
creating_new_EUR_def <- TRUE

PCA_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/09_pca_scores.tsv'
PCA_1KG_SCORES <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_pca_scores_1kg.tsv'

df_PCs <- fread(PCA_SCORES, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df_URVs <- fread(URV_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df_URVs)[names(df_URVs) == 'Samples'] <- 'sample'
df <- merge(df_PCs, df_URVs)
df <- df[sample(nrow(df), replace=FALSE),]

PCs <- c(1,3,5)

if (perform_plotting) {
    for (i in PCs) {
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.PHENOTYPE_COARSE')
        p <- create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1),
          n_x_ticks=5, x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
        p <- p + geom_point(data=df %>% filter(n_URV_SNP > T_nURVSNP | n_URV_indel > T_nURVIndel),
          mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
          inherit.aes=FALSE, shape=4, show.legend=FALSE)
        print(p)
        if(save_figures) {
            ggsave(paste0(PLOTS,'09_PC',i,'_PC',i+1, '_and_URV_marked.jpg'), p, width=144, height=96, units='mm')
        }
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.PROJECT_OR_COHORT')
        create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1, '_batch'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='phenotype.LOCATION')
        create_pretty_scatter(df, aes, save_figure=save_figures, file=paste0(PLOTS,'09_PC',i,'_PC',i+1, '_collection'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}

df_1kg <- fread(PCA_1KG_SCORES, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
    mutate(SUPER_POPULATION = factor(SUPER_POPULATION))
df_1kg <- df_1kg[sample(nrow(df_1kg), replace=FALSE),]
df_1kg <- merge(df_1kg, df_URVs, all.x=TRUE)

if (perform_plotting) {
    for (i in PCs) {
        aes <- aes_string(x=paste0('PC', i), y=paste0('PC', i+1), color='SUPER_POPULATION')
        p <- create_pretty_scatter(df_1kg, aes, #limits=c('Control', 'Bipolar Disorder', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
          add_final_layer=FALSE, #final_layer=subset(df_1kg, Phenotype == 'Bipolar Disorder' | Phenotype=='Control'),
          save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1, '_1kg'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))

        aes <- aes_string(x=paste0('PC', i), y=paste0('PC', i+1), color='PROJECT_OR_COHORT')
        create_pretty_scatter(df_1kg, aes, add_final_layer=FALSE,
          # final_layer=subset(df_1kg, Phenotype == 'Bipolar Disorder' | Phenotype=='Control'),
          save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_1kg_batch'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))

        aes <- aes_string(x=paste0('PC', i), y=paste0('PC', i+1), color='LOCATION')
        p <- create_pretty_scatter(df_1kg, aes, add_final_layer=FALSE,
          # final_layer=subset(df_1kg, Phenotype == 'Bipolar Disorder' | Phenotype=='Control'),
          save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_1kg_collection'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}


df_train = filter(df_1kg, PHENOTYPE_COARSE=='1KG') %>%
  select(c(SUPER_POPULATION, POPULATION, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

df_predict = filter(df_1kg, PHENOTYPE_COARSE!='1KG') %>%
  select(c(s, POPULATION, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,4))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=df_train[PCs_to_use], y=as.factor(as.character(df_train$SUPER_POPULATION)), ntree=10000)
rf_probs <- predict(rf, df_predict[PCs_to_use], type='prob')

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, T_European_RF)
classification <- as.character(predict(rf, df_predict[PCs_to_use]))
df_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
df_predict$classification_strict <- as.factor(classification)

df_classify <- df_predict %>% select(s, classification_strict, classification_loose) %>% inner_join(df_1kg, by='s')

# Include option to use an existing classification

if (!creating_new_EUR_def) {
    # Read in the previously strictly defined Europeans
    european_samples <- fread(EUROPEAN_SAMPLES_STRICT, sep='\t', stringsAsFactors=FALSE, header=FALSE, data.table=FALSE)
    df_classify$classification_strict <- ifelse(df_classify$s %in% european_samples$V1, 'European', 'Non-European')
}

if (perform_plotting) {
    for (i in PCs) {

        if (creating_new_EUR_def)
        {
            aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_loose')
            p <- create_pretty_scatter(df_classify, aes,
              save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_classify_EUR_loose'), n_x_ticks=5,
              x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
            p <- p + geom_point(data=df_1kg %>% filter(SUPER_POPULATION == "EUR"),
              mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
              inherit.aes=FALSE, shape=4, show.legend=FALSE)
            print(p)
        }

        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
        p <- create_pretty_scatter(df_classify, aes,
          save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_classify_EUR_strict'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}

if (creating_new_EUR_def) { 
    df_out_classify_strict <- df_classify %>% filter(classification_strict == 'EUR') %>% select(s)
    df_out_classify_loose <- df_classify %>% filter(classification_loose == 'EUR') %>% select(s)

    # The number of non-European samples
    print(table((df_classify %>% filter(classification_strict != 'EUR'))$phenotype.PHENOTYPE_COARSE))
    print(table((df_classify %>% filter(classification_loose != 'EUR'))$phenotype.PHENOTYPE_COARSE))

    # Print the number of remaining samples.
    print(paste0("Number of European samples using strict classifier: ", dim(df_out_classify_strict)[1]))
    print(paste0("Number of European samples using loose classifier: ", dim(df_out_classify_loose)[1]))

    # Want the number of samples that are non-European and the number of samples that have an excess of URVs.
    fwrite(df_out_classify_strict, file = EUROPEAN_SAMPLES_STRICT, quote=FALSE, row.names=FALSE, col.names=FALSE)
    fwrite(df_out_classify_loose, file = EUROPEAN_SAMPLES_LOOSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Count the number of samples in the strict and non-strict European subset
df_strict <- fread(EUROPEAN_SAMPLES_STRICT, header=FALSE)
df_loose <- fread(EUROPEAN_SAMPLES_LOOSE, header=FALSE)

df_sample_PCs <- data.table(Filter = c("Initial samples",
                          "Non-European, strict filter",
                          "Non-European, loose filter",
                          "Samples after strict European filter"),
                     Samples = c(nrow(df_PCs),
                              nrow(df_PCs) - nrow(df_strict),
                              nrow(df_PCs) - nrow(df_loose),
                              nrow(df_strict)))

fwrite(df_sample_PCs, file='../../samples_BipEx/10_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# Plotting Ultra rare variants against PCs, *before* removal of outliers.

# Plotting the relationship between singleton count and PCs in this dataset.
df$binsPC1 <- cut(df$PC1, breaks = 20)
df$binsPC2 <- cut(df$PC2, breaks = 20)
df$binsPC3 <- cut(df$PC3, breaks = 20)
df$binsPC4 <- cut(df$PC4, breaks = 20)
df$binsPC5 <- cut(df$PC5, breaks = 20)
df$binsPC6 <- cut(df$PC6, breaks = 20)

# Points:
PC1 <- ggplot(df, aes(x=binsPC1, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC2 <- ggplot(df, aes(x=binsPC2, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
PC3 <- ggplot(df, aes(x=binsPC3, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC4 <- ggplot(df, aes(x=binsPC4, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
PC5 <- ggplot(df, aes(x=binsPC5, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC6 <- ggplot(df, aes(x=binsPC6, y=n_URV_SNP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

PC1_Ind <- ggplot(df, aes(x=binsPC1, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC2_Ind <- ggplot(df, aes(x=binsPC2, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
PC3_Ind <- ggplot(df, aes(x=binsPC3, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC4_Ind <- ggplot(df, aes(x=binsPC4, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
PC5_Ind <- ggplot(df, aes(x=binsPC5, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
PC6_Ind <- ggplot(df, aes(x=binsPC6, y=n_URV_indel)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") + coord_flip() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

PC1PC2 <- ggplot(df, aes(PC1, PC2)) + geom_hex(bins=20, show.legend=FALSE) + scale_fill_gradient(name = "count", trans = "log")
PC3PC4 <- ggplot(df, aes(PC3, PC4)) + geom_hex(bins=20, show.legend=FALSE) + scale_fill_gradient(name = "count", trans = "log")
PC5PC6 <- ggplot(df, aes(PC5, PC6)) + geom_hex(bins=20, show.legend=FALSE) + scale_fill_gradient(name = "count", trans = "log")
# Empty plot
empty <- ggplot() + geom_point(aes(1,1), colour="white") +
         theme(axis.ticks=element_blank(), panel.background=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())

grid.arrange(PC1_Ind, empty, empty, PC1, empty, empty, PC1PC2, PC2,  PC2_Ind, ncol=3, nrow=3, widths=c(4, 1, 1), heights=c(1, 1, 4))
grid.arrange(PC3_Ind, empty, empty, PC3, empty, empty, PC3PC4, PC4,  PC4_Ind, ncol=3, nrow=3, widths=c(4, 1, 1), heights=c(1, 1, 4))
grid.arrange(PC5_Ind, empty, empty, PC5, empty, empty, PC5PC6, PC6,  PC6_Ind, ncol=3, nrow=3, widths=c(4, 1, 1), heights=c(1, 1, 4))




