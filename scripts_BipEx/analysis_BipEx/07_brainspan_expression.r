library(data.table)
library(dplyr)

dt_exp <- fread("brainspan_genes_matrix_csv/expression_matrix.csv") %>% rename(row_num = V1)
dt_genes <- fread("brainspan_genes_matrix_csv/rows_metadata.csv") %>% select(row_num, gene_symbol)
dt <- merge(dt_exp, dt_genes, by="row_num")
dt_samples <- fread("brainspan_genes_matrix_csv/columns_metadata.csv")

# Filter to the different ages

get_gene <- function(dt, gene, dt_samples) {
	dt_tmp <- dt %>% filter(gene_symbol == gene) %>% select(-gene_symbol)
	dt <- transpose(data.table(dt_tmp))[-1]
	dt <- cbind(dt, dt_samples) %>% rename(expression = V1)
	dt <- dt %>% mutate(
		age_number = as.integer(sapply(strsplit(age, split=" "), `[`, 1)),
		age_string = sapply(strsplit(age, split=" "), `[`, 2)
		) %>%
		mutate(
			age_bin=case_when(
				(grepl('pcw', age_string) & (age_number <= 12)) ~ 'Early prenatal',
				(grepl('pcw', age_string) & (age_number > 12) & (age_number <= 25)) ~ 'Mid prenatal',
				(grepl('pcw', age_string) & (age_number > 25)) ~ 'Late prenatal',
				(grepl('mos', age_string) | (grepl('yrs', age_string) & (age_number <= 2))) ~ 'Infancy',
				(grepl('yrs', age_string) & (age_number > 2) & (age_number <= 12)) ~ 'Childhood',
				(grepl('yrs', age_string) & (age_number > 12) & (age_number <= 18)) ~ 'Adolescence',
				(grepl('yrs', age_string) & (age_number > 18)) ~ 'Adulthood'
			)
		)

	dt_age <- dt %>% group_by(age_bin) %>% summarise(
		mean=mean(expression),
		median=median(expression),
		lower=quantile(expression, probs=c(0.25)),
		upper=quantile(expression, probs=c(0.75))
		)

	return(dt)
}

dt_plot <- get_gene(dt, 'AKAP11', dt_samples)
p <- ggplot(dt_plot, aes(age_bin, expression)) + 
	geom_boxplot() +
	scale_x_discrete(limits=c(
		'Early prenatal',
		'Mid prenatal',
		'Late prenatal',
		'Infancy',
		'Childhood',
		'Adolescence',
		'Adulthood')
	)
p