library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(logistf)
library(RColorBrewer)
library(pryr)

create_forest_names <- function(dt, BP=TRUE) 
{
	if (BP) {
		dt$Forest_location <- gsub(".*\\, ", "", dt$phenotype.LOCATION)
		dt$Forest_location[grep("Umea, SWE", dt$phenotype.LOCATION)] <- "SWE, Umea"
		dt$Forest_location[grep("Stockholm, SWE", dt$phenotype.LOCATION)] <- "SWE, Stockholm"
		dt$Forest_location[grep("UK|IRE", dt$phenotype.LOCATION)] <- "UK/Ireland"
	} else {
		dt$Forest_location <- gsub(".*\\, ", "", dt$phenotype.LOCATION)
		dt$Forest_location[grep("UK|IRE", dt$phenotype.LOCATION)] <- "UK/Ireland"
	}

	return(dt)
}

tidy_data <- function(dt)
{
	names(dt) <- gsub("pca\\.", "", names(dt))
	names(dt) <- gsub("imputesex\\.impute_sex\\.", "", names(dt))
	return(dt)
}

get_tests_and_covariates <- function(dt)
{
	# Next, want to do both the Firth regression and our regression across all the strata.
	count_tests <- names(dt)[grep("burden", names(dt))]
	count_tests_covariates <- count_tests
	count_tests_covariates <- gsub("_damaging_missense", "", count_tests_covariates)
	count_tests_covariates <- gsub("_other_missense", "", count_tests_covariates)
	count_tests_covariates <- gsub("_non_coding", "", count_tests_covariates)
	count_tests_covariates <- gsub("_coding", "", count_tests_covariates)
	count_tests_covariates <- gsub("_PTV", "", count_tests_covariates)
	count_tests_covariates <- gsub("_synonymous", "", count_tests_covariates)
	count_tests_covariates <- gsub("_MPC_2", "", count_tests_covariates)
	count_tests_covariates <- gsub("_MPC_3", "", count_tests_covariates)
	count_tests_covariates <- gsub("_missense", "", count_tests_covariates)

	where <- c(grep("n_URV$", count_tests), grep("n_URV_SNP$", count_tests), grep("n_URV_indel$", count_tests))
	count_tests <- count_tests[-where]
	count_tests_covariates <- count_tests_covariates[-where]

	where <- c(grep("indel", count_tests), grep("SNP", count_tests))
	count_tests <- count_tests[-where]
	count_tests_covariates <- count_tests_covariates[-where]
	count_tests_covariates_coding <- gsub("n_URV", "n_coding_URV", count_tests_covariates)

	return(list(tests=count_tests, covariates=count_tests_covariates,
		covariates_coding=count_tests_covariates_coding))
}

create_large_forest_plot <- function(
	full_df, out, count_tests, xlabel='', tick_increment=0.1, zero=1,
	clip = c(-Inf, Inf),  pretitle=""
	)
{
	cols <- brewer.pal(n = 8, name = 'Set3')
	cols[2] <- cols[length(cols)]
	cols[length(cols)] <- 'black'
	own <- fpTxtGp(ticks = gpar(cex=1), xlab= gpar(cex=1))

	# Remove the 'coding' URVs these were just an alternate covariate to include
	full_df <- full_df[-grep('n_coding_URV', full_df$variants),]
	# full_df <- full_df[-grep('MPC', full_df$variants),]
	# Now, restrict to the collections of variants and create forest plots.

	# Renaming
	variants <- c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV", "MPC > 2, missense", "MPC > 2, damaging missense")
	ordering <- c("URV_non_coding", "URV_synonymous", "URV_other_missense", "URV_damaging_missense", "URV_PTV", "n_URV_MPC_2_missense", "n_URV_MPC_2_damaging_missense")

	reorder <- function(df, ordering, column) {
		new_ordering <- rep(0, length(ordering))
		k <- 1
		for(i in ordering) {
			new_ordering[k] <- grep(i, (df[,column]))
			k <- k+1
		}
		return(new_ordering)
	}

	titles <- paste0(pretitle, c("Singleton variants", "Singleton variants, not in gnomAD", "Singleton variants, not in gnomAD, pLI > 0.9", "Singleton variants, not in gnomAD, pLI > 0.995"))
	count_tests_matrix <- matrix(count_tests, byrow=TRUE, nrow=4)
	reordering <- reorder(full_df %>% filter(label==names(table(full_df$label))[1]) %>% filter(variants %in% count_tests_matrix[1,]), ordering, 'variants')

	pops <- c("GER", "NED", "SWE, Umea", "USA", "SWE, Stockholm", "UK/Ireland", "Meta-analysis", "All")
	pops_in_plot <- intersect(pops, names(table(full_df$label)))
	cols <- cols[which(pops %in% pops_in_plot)]

	create_to_bind <- function(full_df, to_loop, k)
	{
		df <- full_df %>% filter(label==to_loop[1]) %>% filter(variants %in% count_tests_matrix[k,])
		mean_bind  <- df$mean[reordering]
		lower_bind <- df$lower[reordering]
		upper_bind <- df$upper[reordering]
		
		for(i in 2:length(to_loop)) {
			df <- full_df %>% filter(label==to_loop[i]) %>% filter(variants %in% count_tests_matrix[k,])
			mean_bind  <- cbind(mean_bind,  df$mean[reordering])
			lower_bind <- cbind(lower_bind, df$lower[reordering])
			upper_bind <- cbind(upper_bind, df$upper[reordering])
		}
		box_bind <- 0.01*sqrt(1/(upper_bind-lower_bind))

		return(list(mean_bind=mean_bind, lower_bind=lower_bind, upper_bind=upper_bind, box_bind=box_bind))
	}

	pdf(out, height=8, width=10)
	k <- 1
	for(title in titles) {
		# Define xticks based on range of the data.
		lowest <- min((full_df %>% filter(variants %in% count_tests_matrix[k,]))$lower)
		highest <- max((full_df %>% filter(variants %in% count_tests_matrix[k,]))$upper)
		xlim = c(floor(lowest/tick_increment) * tick_increment, ceiling(highest/tick_increment) * tick_increment)
		for_plot <- create_to_bind(full_df, pops_in_plot, k)

		forestplot(as.matrix(variants),
			for_plot$mean_bind,
			for_plot$lower_bind,
			for_plot$upper_bind,
			lwd_ci=1,
			boxsize = for_plot$box_bind,
			zero=zero,
			xticks=seq(max(clip[1],xlim[1]), min(clip[2], xlim[2]), by=tick_increment),
			xlab=xlabel,
			legend = pops_in_plot,
			col =fpColors(box=cols, line=cols),
			txt_gp=own, title=title, new_page=ifelse(k==1, FALSE, TRUE),
			clip=clip
		)

		k <- k+1
	}

	dev.off()
}

create_small_forest_plot <- function(
	df, variant_column_name, out, xlabel='', tick_increment=0.1, zero=1,
	scale_boxes=1, table_cols=NULL, table_col_names=NULL, clip= c(-Inf, Inf), pretitle="", include_MPC_3=TRUE
	)
{
	cols <- brewer.pal(n = 8, name = 'Set3')
	own <- fpTxtGp(ticks = gpar(cex=1), xlab= gpar(cex=1))
	names(df)[which(names(df) == variant_column_name)] <- 'variants'
	# Remove the 'coding' URVs these were just an alternate covariate to include
	df <- df[-grep('n_coding_URV', df$variants),]

	# Remove the MPC 3 burden - not sufficient samples for a strong signal.
	if (!include_MPC_3) {
		df <- df[-grep('MPC_3', df$variants),]
		
		# Now, restrict to the collections of variants and create forest plots.
		# Renaming
		variants <- c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV", "MPC > 2, missense", "MPC > 2, damaging missense")
		ordering <- c("URV_non_coding", "URV_synonymous", "URV_other_missense", "URV_damaging_missense", "URV_PTV", "URV_MPC_2_missense", "URV_MPC_2_damaging_missense")
	
	} else {
		# Renaming
		variants <- c("Non-coding", "Synonymous", "Other missense", "Damaging missense", "PTV", "MPC > 2, missense", "MPC > 2, damaging missense", "MPC > 3, missense", "MPC > 3, damaging missense")
		ordering <- c("URV_non_coding", "URV_synonymous", "URV_other_missense", "URV_damaging_missense", "URV_PTV", "URV_MPC_2_missense", "URV_MPC_2_damaging_missense", "URV_MPC_3_missense", "URV_MPC_3_damaging_missense")
	}

	reorder <- function(df, ordering, column) {
		new_ordering <- rep(0, length(ordering))
		k <- 1
		for(i in ordering) {
			new_ordering[k] <- grep(i, (df[,column]))
			k <- k+1
		}
		return(new_ordering)
	}

	titles <- paste0(pretitle, c("Singleton variants", "Singleton variants, not in gnomAD", "Singleton variants, not in gnomAD, pLI > 0.9", "Singleton variants, not in gnomAD, pLI > 0.995"))
	count_tests_matrix <- matrix(df$variants, byrow=TRUE, nrow=4)
	reordering <- reorder(df %>% filter(variants %in% count_tests_matrix[1,]), ordering, 'variants')

	pdf(out, height=4, width=10)
	k <- 1

	for(title in titles) {
		# Define xticks based on range of the data.
		lowest <- min((df %>% filter(variants %in% count_tests_matrix[k,]))$lower)
		highest <- max((df %>% filter(variants %in% count_tests_matrix[k,]))$upper)
		xlim = c(floor(lowest/tick_increment) * tick_increment, ceiling(highest/tick_increment) * tick_increment)
		if(!is.null(table_cols)) {
			table_text  <- cbind(as.matrix(variants), (df %>% filter(variants %in% count_tests_matrix[k,]))[reordering, table_cols])
			summary_info <- c(TRUE, rep(FALSE, nrow(table_text)))
			for(col in table_cols) {
				if(is.numeric(table_text[,col]))
					table_text[,col] <- as.character(round(table_text[,col],3))
				table_text[,col] <- as.character(table_text[,col])
			}
			table_text[,1] <- as.character(table_text[,1])
			table_text <- rbind(table_col_names, table_text)
		} else {
			table_text <- as.matrix(variants)
			summary_info <- c(TRUE, rep(FALSE, nrow(table_text)))
			table_text <- rbind("Variants", as.character(table_text))
		}
		forestplot(table_text,
			c(NA, (df %>% filter(variants %in% count_tests_matrix[k,]))$mean[reordering]),
			c(NA, (df %>% filter(variants %in% count_tests_matrix[k,]))$lower[reordering]),
			c(NA, (df %>% filter(variants %in% count_tests_matrix[k,]))$upper[reordering]),
			boxsize = scale_boxes*0.1*sqrt(1/((df %>% filter(variants %in% count_tests_matrix[k,]))$upper[reordering] - 
								   (df %>% filter(variants %in% count_tests_matrix[k,]))$lower[reordering])),
			zero=zero, xticks=seq(xlim[1], xlim[2], by=tick_increment), xlab=xlabel,
			col =fpColors(box=cols, line=cols),
			txt_gp=own, title=title,
			is.summary=summary_info, new_page=ifelse(k==1, FALSE, TRUE),
			clip=clip)
		k <- k+1
	}

	dev.off()
}


get_p <- function(p_val)
{
	p_list_sig <- as.numeric(p_val)
	p_sig <- rep("", length(p_list_sig))
	p_sig[which(p_list_sig < 0.05)] <- "*"
	p_sig[which(p_list_sig < 0.01)] <- "**"
	p_sig[which(p_list_sig < 0.001)] <- "***"
	p_sig[which(p_list_sig < 0.0001)] <- "****"

	return(p_sig)
}

test_burden <- function(
	df, y, x, include_PCs=TRUE, include_sex=TRUE, covariates=NULL,
	poisson=TRUE, logistic=FALSE
	)
{
	formula <- paste0(y, '~', x)
	if(include_PCs) formula <- paste0(formula, '+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10')
	if(include_sex) formula <- paste0(formula, '+ is_female')
	if(!is.null(covariates)) formula <- paste0(formula, '+', paste(covariates, collapse='+'))

	cat('formula evaluating:', formula, '\n')
	formula <- formula(formula)

	if(poisson) {
		fit <- glm(formula=formula, family=poisson(), data=df)
	} else 
	{	
		if(logistic) {
			fit <- glm(formula=formula, family=binomial(), data=df)
		} else {
			fit <- lm(formula=formula, data=df)	
		}
	}

	summ <- summary(fit)
	print(summ)
	x_summ <- ifelse(class(unlist(df[, ..x]))=="logical", paste0(x,"TRUE"), x)
	p_val <- summ$coefficients[x_summ, grep("Pr\\(", colnames(summ$coefficients))]
	variance <- diag(vcov(fit))[x_summ]
	beta = summ$coefficients[x_summ, 'Estimate']
	return(list(fit=fit, Z=sign(beta)*abs(qnorm(p_val/2)), w=sqrt(nrow(df)), beta=beta, var=variance, p_val=p_val))
}

run_Firth_locations_regression <- function(
	dt, count_tests, count_tests_covariates, pheno_var,
	location_var="Forest_location"
	)
{
	dt_forest_Firth <- list()
	Locations <- names(table(dt[[location_var]]))

	kk <- 1
	for(j in count_tests) {
		cat("\nCurrent class of variation:", j, "\n\n")

		if(all(dt[j] == 0)) {
			cat('All 0, move to next test\n')
			kk <- kk+1
			next
		}

		k <- 1
		sum_w_Z <- 0
		sum_w_2 <- 0
		sum_beta_var <- 0
		sum_var <- 0

		# Create an empty data.frame
		dt_forest_Firth[[j]] <- data.table(label=c(Locations, 'Meta-analysis', 'All'))

		for(i in Locations) {
			cat(i, "... ")
			result <- get_Firth_regression(dt, j, pheno_var,
				dataset_location=i, covariates=count_tests_covariates[kk])
			sum_w_Z <- sum_w_Z + result$w * result$Z
			sum_w_2 <- sum_w_2 + result$w^2
			sum_beta_var <- sum_beta_var + (result$beta / result$var)
			sum_var <- sum_var + (1 / result$var)
			dt_forest_Firth[[j]]$mean[k] <- result$beta
			dt_forest_Firth[[j]]$lower[k] <- result$beta -  1.96 * sqrt(result$var)
			dt_forest_Firth[[j]]$upper[k] <- result$beta +  1.96 * sqrt(result$var)
			dt_forest_Firth[[j]]$p_vals[k] <- formatC(result$p_val, digits=3)
			k <- k+1
		}

		# Z-score.
		weighted_Z <- sum_w_Z / sqrt(sum_w_2)
		weighted_beta <- sum_beta_var / sum_var
		se <- sqrt(1 / sum_var)
		# Summary:
		cat('\n\nweighted Z: ', formatC(weighted_Z, digits=3), '\nweighted beta: ',
			formatC(weighted_beta, digits=3), '\nSE: ', formatC(se, digits=3))
		# Associated p-value:
		cat('\np-value: ', formatC(pnorm(weighted_Z, lower.tail=FALSE), digits=3))
		# The meta analysis
		dt_forest_Firth[[j]]$mean[k] <- weighted_beta
		dt_forest_Firth[[j]]$lower[k] <- weighted_beta - 1.96 * se
		dt_forest_Firth[[j]]$upper[k] <- weighted_beta + 1.96 * se
		dt_forest_Firth[[j]]$p_vals[k] <- formatC(2*min(pnorm(weighted_Z, lower.tail=FALSE), pnorm(weighted_Z, lower.tail=TRUE)), digits=3)

		cat("\nFull dataset...")
		result <- get_Firth_regression(dt, j, pheno_var, covariates=count_tests_covariates[kk])

		dt_forest_Firth[[j]]$mean[k+1] <- result$beta
		dt_forest_Firth[[j]]$lower[k+1] <- result$beta -  1.96 * sqrt(result$var)
		dt_forest_Firth[[j]]$upper[k+1] <- result$beta +  1.96 * sqrt(result$var)
		dt_forest_Firth[[j]]$p_vals[k+1] <- formatC(result$p_val, digits=3)

		kk <- kk+1
	}
	return(dt_forest_Firth)
}

get_meta_calcs <- function(fit, meta)
{
	return(
		meta = list(sum_w_Z = meta$sum_w_Z + fit$w * fit$Z,
					sum_w_2 = meta$sum_w_2 + fit$w^2,
					sum_beta_var = meta$sum_beta_var + (fit$beta / fit$var),
					sum_var = meta$sum_var + (1 / fit$var))
	)
}

get_weighted_meta <- function(meta)
{
	# Z-score.
	weighted_Z <- meta$sum_w_Z / sqrt(meta$sum_w_2)
	weighted_beta <- meta$sum_beta_var / meta$sum_var
	se <- sqrt(1 / meta$sum_var)
	# Summary:
	cat('\nweighted Z: ', formatC(weighted_Z, digits=3), '\nweighted beta: ',
		formatC(weighted_beta, digits=3), '\nSE: ', formatC(se, digits=3))
	# Associated p-value:
	cat('\np-value: ', formatC(pnorm(weighted_Z, lower.tail=FALSE), digits=3), '\n')

	return(list(weighted_Z = weighted_Z, weighted_beta=weighted_beta, se=se))
}

run_locations_regression <- function(
	dt, count_tests, pheno_var, count_tests_covariates=NULL,
	location_var="Forest_location"
	) 
{
	dt_forest_logit <- list()
	dt_forest_pois <- list()
	dt_forest_lin <- list()

	Locations <- names(table(dt[[location_var]]))
	dt <- dt %>% mutate(isCase = dt[[pheno_var]])
	
	if(is.null(count_tests_covariates))
		count_tests_covariates <- rep(NULL, length(count_tests))

	kk <- 1
	for(j in count_tests) {
		cat("\nCurrent class of variation:", j, "\n\n")

		if(all(dt[[j]] == 0)) {
			cat('All 0, move to next test\n')
			kk <- kk+1
			next
		}

		# Create an empty data.frame at each entry of the list
		dt_forest_logit[[j]] <- data.table(label=c(Locations, 'Meta-analysis', 'All'))
		dt_forest_pois[[j]] <- data.table(label=c(Locations, 'Meta-analysis', 'All'))
		dt_forest_lin[[j]] <- data.table(label=c(Locations, 'Meta-analysis', 'All'))
		
		k <- 1

		meta_pois <- list(sum_w_Z=0, sum_w_2=0, sum_beta_var=0, sum_var=0)
		meta_logit <- list(sum_w_Z=0, sum_w_2=0, sum_beta_var=0, sum_var=0)
		meta_lin <- list(sum_w_Z=0, sum_w_2=0, sum_beta_var=0, sum_var=0)

		for(i in Locations)
		{
			dt_filtered <- dt %>% filter(dt[[location_var]] == i)
			cat(paste0(i, "... "))
			cat(sum(dt_filtered[[j]] > 0), '...')
			if (sum(dt_filtered[[j]] > 0) < 15) {
				cat('Not sufficient data to fit...\n')
				dt_forest_pois[[j]]$mean[k] <- NA
				dt_forest_pois[[j]]$lower[k] <- NA
				dt_forest_pois[[j]]$upper[k] <- NA
				dt_forest_pois[[j]]$p_vals[k] <- NA

				dt_forest_logit[[j]]$mean[k] <- NA
				dt_forest_logit[[j]]$lower[k] <- NA
				dt_forest_logit[[j]]$upper[k] <- NA
				dt_forest_logit[[j]]$p_vals[k] <- NA

				dt_forest_lin[[j]]$mean[k] <- NA
				dt_forest_lin[[j]]$lower[k] <- NA
				dt_forest_lin[[j]]$upper[k] <- NA
				dt_forest_lin[[j]]$p_vals[k] <- NA

				k <- k+1

				next
			}

			# Poisson regression
			fit <- test_burden(dt_filtered, j, 'isCase', covariates=count_tests_covariates[kk])
			meta_pois <- get_meta_calcs(fit, meta_pois)

			conf_int <- suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95))
			dt_forest_pois[[j]]$mean[k] <- fit$beta
			dt_forest_pois[[j]]$lower[k] <- conf_int[1]
			dt_forest_pois[[j]]$upper[k] <- conf_int[2]
			dt_forest_pois[[j]]$p_vals[k] <- formatC(signif(fit$p_val, digits=3))

			# Logistic regression
			fit <- test_burden(dt_filtered, 'isCase', j, logistic=TRUE, poisson=FALSE, covariates=count_tests_covariates[kk])
			meta_logit <- get_meta_calcs(fit, meta_logit)

			conf_int <- suppressMessages(confint(fit$fit, j, level=0.95))
			dt_forest_logit[[j]]$mean[k] <- exp(fit$beta)
			dt_forest_logit[[j]]$lower[k] <- exp(conf_int[1])
			dt_forest_logit[[j]]$upper[k] <- exp(conf_int[2])
			dt_forest_logit[[j]]$p_vals[k] <- formatC(signif(fit$p_val, digits=3))

			# Linear regression
			fit <- test_burden(dt_filtered, j, 'isCase', poisson=FALSE, logistic=FALSE, covariates=count_tests_covariates[kk])
			meta_lin <- get_meta_calcs(fit, meta_lin)

			conf_int <- suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95))
			dt_forest_lin[[j]]$mean[k] <- fit$beta
			dt_forest_lin[[j]]$lower[k] <- conf_int[1]
			dt_forest_lin[[j]]$upper[k] <- conf_int[2]
			dt_forest_lin[[j]]$p_vals[k] <- formatC(signif(fit$p_val, digits=3))

			k <- k+1
		}

		cat("\n\nNow perform Meta-analysis...\n")

		weighted_pois <- get_weighted_meta(meta_pois)
		weighted_logit <- get_weighted_meta(meta_logit)
		weighted_lin <- get_weighted_meta(meta_lin)

		dt_forest_pois[[j]]$mean[k] <- weighted_pois$weighted_beta
		dt_forest_pois[[j]]$lower[k] <- weighted_pois$weighted_beta - 1.96 * weighted_pois$se
		dt_forest_pois[[j]]$upper[k] <- weighted_pois$weighted_beta + 1.96 * weighted_pois$se
		dt_forest_pois[[j]]$p_vals[k] <- formatC(2*min(pnorm(weighted_pois$weighted_Z, lower.tail=FALSE), pnorm(weighted_pois$weighted_Z, lower.tail=TRUE)), digits=3)

		dt_forest_logit[[j]]$mean[k] <- exp(weighted_logit$weighted_beta)
		dt_forest_logit[[j]]$lower[k] <- exp(weighted_logit$weighted_beta - 1.96 * weighted_logit$se)
		dt_forest_logit[[j]]$upper[k] <- exp(weighted_logit$weighted_beta + 1.96 * weighted_logit$se)
		dt_forest_logit[[j]]$p_vals[k] <- formatC(2*min(pnorm(weighted_logit$weighted_Z, lower.tail=FALSE), pnorm(weighted_logit$weighted_Z, lower.tail=TRUE)), digits=3)

		dt_forest_lin[[j]]$mean[k] <- weighted_lin$weighted_beta
		dt_forest_lin[[j]]$lower[k] <- weighted_lin$weighted_beta - 1.96 * weighted_lin$se
		dt_forest_lin[[j]]$upper[k] <- weighted_lin$weighted_beta + 1.96 * weighted_lin$se
		dt_forest_lin[[j]]$p_vals[k] <- formatC(2*min(pnorm(weighted_lin$weighted_Z, lower.tail=FALSE), pnorm(weighted_lin$weighted_Z, lower.tail=TRUE)), digits=3)

		# Poisson regression
		cat("\nFull dataset...")
		fit <- test_burden(dt, j, 'isCase', covariates=count_tests_covariates[kk])
		conf_int <- suppressMessages(suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95)))
		dt_forest_pois[[j]]$mean[k+1] <- fit$beta
		dt_forest_pois[[j]]$lower[k+1] <- conf_int[1]
		dt_forest_pois[[j]]$upper[k+1] <- conf_int[2]
		dt_forest_pois[[j]]$p_vals[k+1] <- formatC(signif(fit$p_val, digits=3))

		# Logistic regression
		fit <- test_burden(dt, 'isCase', j, logistic=TRUE, poisson=FALSE, covariates=count_tests_covariates[kk])
		conf_int <- suppressMessages(confint(fit$fit, j, level=0.95))
		dt_forest_logit[[j]]$mean[k+1] <- exp(fit$beta)
		dt_forest_logit[[j]]$lower[k+1] <- exp(conf_int[1])
		dt_forest_logit[[j]]$upper[k+1] <- exp(conf_int[2])
		dt_forest_logit[[j]]$p_vals[k+1] <- formatC(signif(fit$p_val, digits=3))

		# Linear regression
		fit <- test_burden(dt, j, 'isCase', poisson=FALSE, logistic=FALSE, covariates=count_tests_covariates[kk])
		conf_int <- suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95))
		dt_forest_lin[[j]]$mean[k+1] <- fit$beta
		dt_forest_lin[[j]]$lower[k+1] <- conf_int[1]
		dt_forest_lin[[j]]$upper[k+1] <- conf_int[2]
		dt_forest_lin[[j]]$p_vals[k+1] <- formatC(signif(fit$p_val, digits=3))

		kk <- kk+1
	}

	return(locations_regression = list(logit=dt_forest_logit, lin=dt_forest_lin, pois=dt_forest_pois))
}

run_collection_burden_regression <- function(
	count_tests, dt, pheno_var, count_tests_covariates=NULL, run_firth=FALSE,
	run_linear=TRUE, run_logistic=TRUE, run_poisson=TRUE
	)
{
	dt_forest_pois <- data.frame(label=count_tests)
	dt_forest_logit <- data.frame(label=count_tests)
	dt_forest_lin <- data.frame(label=count_tests)
	dt_forest_firth <- data.frame(label=count_tests)

	if(is.null(count_tests_covariates))
		count_tests_covariates <- rep(NULL, length(count_tests))

	dt <- dt %>% mutate(isCase = dt[[pheno_var]])

	k <- 1
	for(i in count_tests) {
		if(all(dt[, ..i] == 0)) {
			k <- k+1
			next
		}

		# Poisson regression
		print(count_tests_covariates[[k]])
		if (run_poisson) {
			fit <- test_burden(dt, i, 'isCase', covariates=count_tests_covariates[[k]])
			conf_int <- suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95))
			dt_forest_pois$mean[k] <- fit$beta
			dt_forest_pois$lower[k] <- conf_int[1]
			dt_forest_pois$upper[k] <- conf_int[2]
			dt_forest_pois$p_vals[k] <- formatC(signif(fit$p_val, digits=3))
		}

		# Logistic regression
		if (run_logistic) {
			fit <- test_burden(dt, 'isCase', i, poisson=FALSE, logistic=TRUE, covariates=count_tests_covariates[[k]])
			conf_int <- suppressMessages(confint(fit$fit, i, level=0.95))
			dt_forest_logit$mean[k] <- exp(fit$beta)
			dt_forest_logit$lower[k] <- exp(conf_int[1])
			dt_forest_logit$upper[k] <- exp(conf_int[2])
			dt_forest_logit$p_vals[k] <- formatC(signif(fit$p_val, digits=3))
		}

		# Linear regression
		if (run_linear) {
			fit <- test_burden(dt, i, 'isCase', poisson=FALSE, logistic=FALSE, covariates=count_tests_covariates[[k]])
			conf_int <- suppressMessages(confint(fit$fit, paste0('isCase', 'TRUE'), level=0.95))
			dt_forest_lin$mean[k] <- fit$beta
			dt_forest_lin$lower[k] <- conf_int[1]
			dt_forest_lin$upper[k] <- conf_int[2]
			dt_forest_lin$p_vals[k] <- formatC(signif(fit$p_val, digits=3))
		}

		# Firth regression
		if (run_firth) {
			fit <- get_Firth_regression(dt, i, pheno_var, covariates=count_tests_covariates[[k]])
			dt_forest_firth$mean[k] <- fit$beta
			dt_forest_firth$lower[k] <- fit$beta -  1.96 * sqrt(fit$var)
			dt_forest_firth$upper[k] <- fit$beta +  1.96 * sqrt(fit$var)
			dt_forest_firth$p_vals[k] <- formatC(fit$p_val, digits=3)
		}

		k <- k+1
	}

	list_to_return <- list()

	if (run_firth) list_to_return$firth <- dt_forest_firth
	if (run_poisson) list_to_return$pois <- dt_forest_pois
	if (run_logistic) list_to_return$logit <- dt_forest_logit
	if (run_linear) list_to_return$lin <- dt_forest_lin

	return(list_to_return)

}

naive_loop <- function(dt)
{
	genes <- names(table(dt$gene_symbol))
	gene <- "ST3GAL2"
	phenotype <- "is_BPSCZ"
	consequence_categories <- names(table(dt$consequence_category))
	cols_to_extract <- paste0(phenotype, c(".case_count", ".case_no_count", ".control_count", ".control_no_count"))
	locations <- names(table(dt$forest_location))
	dt <- dt %>% filter(gene_symbol == gene) %>% select(forest_location, consequence_category, cols_to_extract)

	for(consequence in consequence_categories)
	{
		cat(consequence, "\n")
		dt_tmp <- dt %>% filter(consequence_category == consequence)
		locations <- names(table((dt_tmp)$forest_location))
		
		CMH_matrix <- array(dim=c(2,2,length(locations)))
		k <- 1
		if(length(locations) > 0) {
			for(i in locations) {
				CMH_matrix[,,k] <- matrix(unlist(dt_tmp %>% filter(forest_location == i) %>% select(cols_to_extract)), byrow=TRUE, ncol=2)
				k <- k+1
			}
			# print(CMH_matrix)
			print(mantelhaen.test(CMH_matrix))
		}
		print(fisher.test(apply(CMH_matrix, c(1,2), sum)))
		# print(apply(CMH_matrix, c(1,2), sum))
	}

}

OR_naive <- function(CMH_matrix) {
	numer <- 0
	denom <- 0
	numer_2 <- 0
	denom_2 <- 0
	for(i in 1:dim(CMH_matrix)[3]) {
		
		CMH <- CMH_matrix[,,i]
		T <- sum(CMH)

		A <- CMH[1,1]
		B <- CMH[1,2]
		C <- CMH[2,1]
		D <- CMH[2,2]

		n1 <- A+C
		n2 <- B+D

		m1 <- A+B
		m2 <- C+D

		numer <- numer + (A * D) / T
		denom <- denom + (B * C) / T
		numer_2 <- numer_2 + A - n1*m1/T
		denom_2 <- denom_2 + (n1/T)*(n2/T)*m1*(m2/(T-1)) 
	}

}

get_CMH_genes <- function(phenotype, dt,
	counts_vec=c(".case_count", ".case_no_count", ".control_count", ".control_no_count")
) {
	# contingency table is:
	# | x, B | k (= number of balls drawn (= number of cases))
	# | C, D | l (= number of balls not drawn (= number of controls))
	#   m  n   T (= total balls in the urn (= sample size)) 
	# x (notation in rhyper R docs) is number of white balls drawn, 
	# B is number of black balls drawn.
	# C and D are the numbers of white and black balls not drawn.
	# m and n are the numbers of white and black balls in the urn.
	dt <- dt %>% 
	mutate(
		x=dt[[paste0(phenotype, counts_vec[1])]], 
		B=dt[[paste0(phenotype, counts_vec[2])]],
		C=dt[[paste0(phenotype, counts_vec[3])]],
		D=dt[[paste0(phenotype, counts_vec[4])]]
	) %>% 
	filter(!is.na(x)) %>%
	mutate(
		T=x+B+C+D,
		k=x+B,
		l=C+D,
		m=x+C,
		n=B+D
	) %>% 
	mutate(
		mk_T=m*(k/T),
		nmlk_T2T=(n/T)*(m/T)*(l/(T-1))*k
	) %>%
	select(gene_symbol, consequence_category, forest_location, x, B, C, D, T, mk_T, nmlk_T2T) %>% 
	mutate(
		e_numer=x-mk_T,
		R_numer=(x*D)/T,
		R_denom=(B*C)/T
	)

	dt_summary <- dt %>%
	dplyr::group_by(gene_symbol, consequence_category) %>% 
	summarize(
		mean_R_numer=mean(R_numer),
		mean_R_denom=mean(R_denom),
		sum_e_numer=sum(e_numer),
		sum_e_denom=sum(nmlk_T2T)
	)

	dt_summary <- data.table(dt_summary) %>% 
	mutate(R=mean_R_numer/mean_R_denom) %>% 
	mutate(
		# No continuity correction
		e_CMH=ifelse(is.nan(R), NA, (sum_e_numer^2)/sum_e_denom),
		# Continuity correction
		e_CMH_cc=ifelse(is.nan(R), NA, (abs(sum_e_numer)-0.5)^2/sum_e_denom)
	)
	dt_summary <- data.table(dt_summary)
	dt_summary <- dt_summary[,pval:=pchisq(e_CMH, 1, lower=FALSE)]
	dt_summary <- dt_summary[,pval_cc:=pchisq(e_CMH_cc, 1, lower=FALSE)]

	dt_pval <- dt_summary %>% select(gene_symbol, consequence_category, pval, pval_cc, R, e_CMH_cc, e_CMH) %>% arrange(consequence_category, pval_cc) 
	return(dt_pval)
}

create_permutation_CMH_genes <- function(
	n_permutations, phenotype, dt, counts_vec=c(".case_count", ".case_no_count", ".control_count", ".control_no_count")
) {

	# contingency table is:
	# | x, B | k (= number of balls drawn (= number of cases))
	# | C, D | l (= number of balls not drawn (= number of controls))
	#   m  n   T (= total balls in the urn (= sample size)) 
	# x (notation in rhyper R docs) is number of white balls drawn, 
	# B is number of black balls drawn.
	# C and D are the numbers of white and black balls not drawn.
	# m and n are the numbers of white and black balls in the urn.
	dt <- dt %>% 
		mutate(
			x=dt[[paste0(phenotype, counts_vec[1])]], 
			B=dt[[paste0(phenotype, counts_vec[2])]],
			C=dt[[paste0(phenotype, counts_vec[3])]],
			D=dt[[paste0(phenotype, counts_vec[4])]]) %>% 
		filter(!is.na(x)) %>%
		mutate(
			T=x+B+C+D,
			k=x+B,
			l=C+D,
			m=x+C,
			n=B+D) %>% 
		mutate(
			mk_T=m*(k/T),
			nmlk_T2T=(n/T)*(m/T)*(l/(T-1))*k) %>% 
		select(gene_symbol, consequence_category, forest_location, n, m, l, k, T, mk_T, nmlk_T2T)

	dt <- data.table(dt)

	sample_contingency <- function(m, n, k, T, mk_T)
	{
		x <- mapply(rhyper, 1, m, n, k)
		B <- k-x
		C <- m-x
		D <- n-B

		R_numer <- (x*D)/T
		R_denom <- (B*C)/T

		e_numer <- x - mk_T

		return(list(e_numer=e_numer, R_numer=R_numer, R_denom=R_denom))
	}

	init <- TRUE
	for(i in 1:n_permutations)
	{
		dt[, c("e_numer", "R_numer", "R_denom") := sample_contingency(m, n, k, T, mk_T)]
		# quartz()
		# hist(log(dt$R_denom))
		dt_summary <- dt %>% 
		group_by(gene_symbol, consequence_category) %>% 
		summarize(
			mean_R_numer=mean(R_numer),
			mean_R_denom=mean(R_denom),
			sum_e_numer=sum(e_numer),
			sum_e_denom=sum(nmlk_T2T)
		)
		dt_summary <- data.table(dt_summary) %>% 
		mutate(R=mean_R_numer/mean_R_denom) %>% 
		mutate(
			# No continuity correction
			e_CMH=ifelse(is.nan(R), NA, (sum_e_numer^2)/sum_e_denom),
			# Continuity correction
			e_CMH_cc=ifelse(is.nan(R), NA, (abs(sum_e_numer)-0.5)^2/sum_e_denom)
		)

		dt_summary <- data.table(dt_summary)
		dt_summary <- dt_summary[,pval:=pchisq(e_CMH, 1, lower=FALSE)]
		dt_summary <- dt_summary[,pval_cc:=pchisq(e_CMH_cc, 1, lower=FALSE)]

		cat('permutation =', i, '\n')
		if (init) {
			dt_pval <- dt_summary %>% arrange(consequence_category, pval_cc) %>% select(gene_symbol, consequence_category, pval_cc, R)
			init <- FALSE
		} else {
			dt_tmp <- dt_summary %>% arrange(consequence_category, pval_cc) %>% select(gene_symbol, consequence_category, pval_cc, R)
			dt_pval <- dt_pval %>% mutate(R = R + dt_tmp[['R']], pval_cc = pval_cc + dt_tmp[['pval_cc']])
		}
	}

	dt_pval <- data.table(dt_pval)
	dt_pval[, pval_cc:=(pval_cc/n_permutations)]
	dt_pval[, R:=R/n_permutations]

	return(dt_pval)
}


get_fisher_genes <- function(phenotype, dt,
	counts_vec=c(".case_count", ".case_no_count", ".control_count", ".control_no_count")
) {
	# contingency table is:
	# | x, B | k (= number of balls drawn (= number of cases))
	# | C, D | l (= number of balls not drawn (= number of controls))
	#   m  n   T (= total balls in the urn (= sample size)) 
	# x (notation in rhyper R docs) is number of white balls drawn, 
	# B is number of black balls drawn.
	# C and D are the numbers of white and black balls not drawn.
	# m and n are the numbers of white and black balls in the urn.
	dt <- dt %>% 
	mutate(
		x=dt[[paste0(phenotype, counts_vec[1])]], 
		B=dt[[paste0(phenotype, counts_vec[2])]],
		C=dt[[paste0(phenotype, counts_vec[3])]],
		D=dt[[paste0(phenotype, counts_vec[4])]]) %>% 
	select(gene_symbol, consequence_category, x, B, C, D) %>% 
	filter(!is.na(x))

	fisher <- function(x, B, C, D)
	{
		result <- fisher.test(matrix(c(x, C, B, D), nrow=2, byrow=FALSE))
		return(list(pval=result$p.value, OR=result$estimate))
	}

	list_fisher <- function(x, B, C, D) {
		result <- mapply(fisher, x, B, C, D)
		return(list(pval=unlist(result[1,]), OR=unlist(result[2,])))
	}

	dt <- data.table(dt)
	dt[, c("pval","OR"):=list_fisher(x, B, C, D)]
	dt_pval <- dt %>% select(gene_symbol, consequence_category, pval, OR) %>% arrange(consequence_category, pval)

	return(dt_pval)
}


create_permutation_genes <- function(
	n_permutations, phenotype, dt, counts_vec=c(".case_count", ".case_no_count", ".control_count", ".control_no_count")
) {

	# contingency table is:
	# | x, B | k (= number of balls drawn (= number of cases))
	# | C, D | l (= number of balls not drawn (= number of controls))
	#   m  n   T (= total balls in the urn (= sample size)) 
	# x (notation in rhyper R docs) is number of white balls drawn, 
	# B is number of black balls drawn.
	# C and D are the numbers of white and black balls not drawn.
	# m and n are the numbers of white and black balls in the urn.

	dt <- dt %>% 
		mutate(
			x=dt[[paste0(phenotype, counts_vec[1])]], 
			B=dt[[paste0(phenotype, counts_vec[2])]],
			C=dt[[paste0(phenotype, counts_vec[3])]],
			D=dt[[paste0(phenotype, counts_vec[4])]]) %>% 
		mutate(
			k=x+B,
			m=x+C,
			n=B+D) %>% filter(!is.na(x)) %>%
		select(gene_symbol, consequence_category, n, m, k)

	dt <- data.table(dt)

	sample_contingency <- function(m, n, k)
	{
		x <- mapply(rhyper, 1, m, n, k)
		B <- k-x
		C <- m-x
		D <- n-B

		fisher <- function(x, B, C, D) {
			return(fisher.test(matrix(c(x, C, B, D), nrow=2, byrow=FALSE))$p.value)
		}

		pval <- mapply(fisher, x, B, C, D)

		return(pval)
	}

	init <- TRUE
	for(i in 1:n_permutations)
	{
		dt[, pval:=sample_contingency(m, n, k)]

		cat('permutation =', i, '\n')
		if (init) {
			dt_pval <- dt %>% arrange(consequence_category, pval) %>% select(gene_symbol, consequence_category, pval)
			init <- FALSE
		} else {
			dt_tmp <- dt %>% arrange(consequence_category, pval) %>% select(gene_symbol, consequence_category, pval)
			dt_pval <- dt_pval %>% mutate(pval = pval + dt_tmp[['pval']])
		}
	}

	dt_pval <- data.table(dt_pval)
	dt_pval[, pval:=(pval/n_permutations)]

	return(dt_pval)
}

get_contingency <- function(
	dt, variant_collection, pheno_var, location_var = "Forest_location",
	dataset_location=NULL
	)
{
	if(is.null(dataset_location)) {
		dt_cont <- dt %>% mutate(isCase = dt[[pheno_var]],
			varCollection = dt[[variant_collection]])
	} else {
		dt_cont <- dt %>% filter(dt[[location_var]] == dataset_location)
		dt_cont <- dt_cont %>% mutate(isCase = dt_cont[[pheno_var]],
			varCollection = dt_cont[[variant_collection]])
	}

	dt_cont <- dt_cont %>% mutate(carrier = (varCollection > 0)) %>% select(carrier, isCase)
	cont_table <- matrix(
		c(sum((dt_cont %>% filter(isCase))$carrier),
		  sum(!(dt_cont %>% filter(isCase))$carrier),
		  sum((dt_cont %>% filter(!isCase))$carrier),
		  sum(!(dt_cont %>% filter(!isCase))$carrier)),
		nrow=2, ncol=2, byrow=TRUE)

	# Sanity check that these sum to the total number of samples in that location.
	return(cont_table)
}

run_collection_CMH <- function(
	dt, count_tests, pheno_var, location_var = "Forest_location", 
	verbose=FALSE, just_Fisher=FALSE, exact=TRUE
	)
{
	Locations <- names(table(dt[[location_var]]))
	dt_CMH <- data.table(labels=count_tests)
	kk <- 1
	CMH_list <- list()

	for(j in count_tests)
	{
		cat("\n", j, "\n")
		if(!just_Fisher) {
			CMH_matrix <- array(dim=c(2,2,length(Locations)))
			k <- 1

			for(i in Locations) {
				cat(i, "\n")
				CMH_matrix_sub_array <- get_contingency(dt, j, pheno_var, dataset_location=i)

				if(verbose) 
					print(CMH_matrix_sub_array)
				CMH_matrix[,,k] <- CMH_matrix_sub_array
				k <- k+1
			}

			CMH_list[[j]] <- CMH_matrix
			CMH <- mantelhaen.test(CMH_matrix, exact=exact)

			dt_CMH$CMH_p[kk] <- CMH$p.value
			dt_CMH$CMH_OR[kk] <- CMH$estimate
		}

		cat("All locations\n")
		cont_all <- get_contingency(dt, j, pheno_var)
		Fish <- fisher.test(cont_all)
		dt_CMH$Fisher_p[kk] <- Fish$p.value
		dt_CMH$Fisher_OR[kk] <- Fish$estimate
		kk <- kk+1
	}
	return(dt_CMH)
}

get_Firth_regression <- function(
	df, variant_collection, pheno_var, location_var = "Forest_location",
	dataset_location=NULL, include_PCs=TRUE, include_sex=TRUE, covariates=NULL
	) 
{

	if(is.null(dataset_location)) {
		df <- df %>% mutate(isCase = df[[pheno_var]])
	} else {
		df <- df %>% filter(df[[location_var]] == dataset_location)
		df <- df %>% mutate(isCase = df[[pheno_var]])
	}

	formula <- paste0('isCase', '~', variant_collection)

	if(include_PCs) formula <- paste0(formula, '+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10')
	if(include_sex) formula <- paste0(formula, '+ is_female')
	if(!is.null(covariates)) formula <- paste0(formula, '+', covariates)

	formula <- formula(formula)

	fit <- logistf(formula=formula, data=df)
	summary <- summary(fit)
	where <- which(names(summary$prob) == variant_collection)
	p_val <- summary$prob[variant_collection]
	variance <- diag(vcov(fit))[where]
	beta = summary$coefficients[variant_collection]

	return(list(Z=sign(beta)*abs(qnorm(p_val/2)), w=sqrt(nrow(df)), beta=beta, var=variance, p_val=p_val))
}

create_fisher <- function(vector_variants_cases_controls, n_cases, n_controls) {
	results <- fisher.test(cbind(vector_variants_cases_controls, (c(n_cases, n_controls) - vector_variants_cases_controls)))
	return(c(pval=results$p.value, lower_CI=results$conf.int[1], upper_CI=results$conf.int[2], OR=results$estimate))
}

create_fisher_lookup_table <- function(
	n_cases, n_controls, singleton_df, case_col, control_col
	)
{	
	maximum_possible_carriers <- max(rowSums(singleton_df %>% select(case_col, control_col)))
	# Predefine contingency table
	contingency_table <- matrix(0, nrow=2, ncol=2)
	Fisher_precompute <- matrix(0,
		nrow=(min(n_cases, maximum_possible_carriers)+1),
		ncol=(min(n_controls, maximum_possible_carriers)+1))
	
	for (i in 1:nrow(Fisher_precompute))
	{
		contingency_table[1,1] <- i - 1
		contingency_table[1,2] <- n_cases - contingency_table[1,1]
		
		for (j in 1:ncol(Fisher_precompute))
		{
			contingency_table[2,1] <- j - 1
			contingency_table[2,2] <- n_controls - contingency_table[2,1]
			Fisher_precompute[i,j] <- fisher.test(contingency_table)$p.value
		}
		cat(paste0(i,' of ', dim(Fisher_precompute)[1],' possibilities \n'))
	}
	return(Fisher_precompute)
}

create_permutation_p_fisher <- function(
	n_permutations, singleton_df, n_cases, n_controls,
	Fisher_lookup, case_col, control_col)
{
	n_genes <- nrow(singleton_df)
	n <- n_cases + n_controls
	p_values <- rep(0, n_genes)
	# For each gene, find out how many mutations need to be allocated.
	mutations_for_allocation <- rowSums(singleton_df %>% select(case_col, control_col))

	apply_rhyper <- function(mutations_for_allocation, n_cases, n_controls, n_perm) {
		rhyper(n_perm, n_cases, n_controls, mutations_for_allocation)
	}

	for(i in 1:n_permutations) {
		n_cases_with_mut <- sapply(mutations_for_allocation, apply_rhyper, n_cases, n_controls, 1)
		n_controls_with_mut <- mutations_for_allocation - n_cases_with_mut
		p_values <- p_values + sort(Fisher_lookup[cbind(n_cases_with_mut+1, n_controls_with_mut+1)])
		if ((i %% 10000) == 0) {
			cat('permutation =',i,'\n')
		}
	}

	p_values <- p_values / n_permutations
	return(p_values)
}

create_permutation_p_fisher <- function(
	n_permutations, dt, n_cases, n_controls,
	Fisher_lookup, case_col, control_col, p_value_obs)
{
	n_genes <- nrow(dt)
	n <- n_cases + n_controls

	# For each gene, find out how many mutations need to be allocated.
	mutations_for_allocation <- rowSums(dt %>% select(case_col, control_col))

	apply_rhyper <- function(mutations_for_allocation, n_cases, n_controls, n_perm) {
		rhyper(n_perm, n_cases, n_controls, mutations_for_allocation)
	}

	dt$p_perm <- 0

	for(i in 1:n_permutations) {
		n_cases_with_mut <- sapply(mutations_for_allocation, apply_rhyper, n_cases, n_controls, 1)
		n_controls_with_mut <- mutations_for_allocation - n_cases_with_mut
		dt$p_perm <- dt$p_perm + (Fisher_lookup[cbind(n_cases_with_mut+1, n_controls_with_mut+1)] <= dt$pval)
		if ((i %% 10000) == 0) {
			cat('permutation =',i,'\n')
		}
	}

	dt <- dt %>% mutate(p_perm = p_perm/n_permutations)
	return(dt)
}

run_qq_create_fisher_and_perm <- function(dt, n_cases, n_controls, 
	case_col, control_col, n_permutations=1000, consequence="ptv",
	save_lookup=TRUE, lookup_filename='lookup.Rdata',
	read_lookup=FALSE, include_col=FALSE)
{	
	cat(paste0("\nn cases: ", n_cases, "\nn controls: ", n_controls, "\n"))
	print(head(dt %>% select(case_col, control_col)))

	dt <- dt %>% filter(consequence_category==consequence)
	results <- t(apply(dt %>% select(case_col, control_col), 1, create_fisher, n_cases, n_controls))
	results_dt <- cbind(dt, results)
	# return(results_dt)

	if (read_lookup) {
		cat("Reading Fisher look-up...\n")
		if (file.exists(lookup_filename)) {
			load(lookup_filename)
		} else {
			cat("Specified Rdata file does not exist, creating Fisher look-up...\n")
			Fisher_lookup <- create_fisher_lookup_table(
				n_cases, n_controls, dt %>% filter(consequence_category==consequence),
				case_col, control_col
			)
		}
	} else {
		cat("Creating Fisher look-up...\n")
		Fisher_lookup <- create_fisher_lookup_table(
			n_cases, n_controls, dt %>% filter(consequence_category==consequence),
			case_col, control_col
		)
	}

	if(save_lookup) {
		cat("Saving Fisher look-up...\n")
		save(Fisher_lookup, file=lookup_filename)
	}

	cat("Creating permutation p-values...\n")
	p_perm <- create_permutation_p_fisher(
		n_permutations, dt %>% filter(consequence_category==consequence),
		n_cases, n_controls, Fisher_lookup, case_col, control_col)
	p_obs <- -log10(results_dt %>% filter(consequence_category == consequence) %>% arrange(pval) %>% select(pval))$pval

	dt_plot <- data.table(
		p_perm = -log10(p_perm),
		pval = p_obs,
		labels = (results_dt %>% filter(consequence_category==consequence) %>%
			arrange(pval) %>% select(gene_symbol))$gene_symbol)

	return(list(dt_plot=dt_plot, results_dt=results_dt))
}

create_qq_plot_and_table <- function(dt_plot_filename, qq_labels=20, n_out=2000,
	plot_title="", table_out="table", include_col=TRUE) 
{
	load(dt_plot_filename)
	p_obs <- dt_plot$pval
	dt_plot <- merge(dt_plot, results_dt %>% rename(labels=gene_symbol) %>% select(-pval))
	dt_plot <- dt_plot %>% rename(odds_ratio = `OR.odds ratio`) %>% 
		mutate(log_OR = log10(odds_ratio))
    max_replace <- max(abs(dt_plot$log_OR[is.finite(dt_plot$log_OR)]))
    dt_plot$log_OR[!is.finite(dt_plot$log_OR)] <- sign(dt_plot$log_OR[!is.finite(dt_plot$log_OR)]) * max_replace
	cat("Creating QQ plot...\n")
	if (include_col) {
		p <- create_pretty_qq_plot(
			dt_plot, aes(x=p_perm, y=pval, color=log_OR),
			n_to_include=qq_labels, cex_label=2,
			plot_title=plot_title, gradient=include_col, 
			gradient_title=TeX("$\\log(OR)$"))
	} else {
		p <- create_pretty_qq_plot(
			dt_plot, aes(x=p_perm, y=pval),
			n_to_include=qq_labels, cex_label=2,
			plot_title=plot_title, gradient=include_col)
	}
	cat("Created QQ plot...\n")
	cat("Restricting to the top", n_out, "p-values...\n")
	results_out_dt <- (results_dt  %>% arrange(pval))[1:n_out,]
	cat("Writing to .tsv file...\n")
	fwrite(results_out_dt[1:min(nrow(results_out_dt), n_out),], file=paste0(table_out, '.tsv'))
	cat("Done.\n\n")
	return(p)
}

create_burden_file <- function(dt, gene_list, consequence_categories, file_out_prefix="BPSCZ")
{
	gene_list <- gene_list %>% rename(gene_symbol = gene)
	gene_sets <- unique(gene_list$geneset_name)
	gene_list <- data.table(gene_list)
	setkey(gene_list, 'gene_symbol')

	for(consequence in consequence_categories)
	{
		cat("Consequence_category:", consequence, "\n")
		dt_cons <- dt %>% filter(consequence_category==consequence)
		dt_cons <- dt_cons %>% filter(apply(dt_cons[, -c(1,2)], 1, function(x){any(x) > 0}))
		dt_cons <- data.table(dt_cons)
		dt_cons[,consequence_category:=NULL]
		setkey(dt_cons, 'gene_symbol')

		n_per_loop <- 400
		n_loops <- ceiling((ncol(dt_cons)-1)/n_per_loop)
		dt_list <- list()
		start <- 2
		for(i in 1:n_loops) {
			cat("Current loop:", i, "\n")
			cols <- names(dt_cons)[c(1, seq(start, min((start+n_per_loop), ncol(dt_cons))))]
			dt_tmp <- dt_cons[, ..cols]
			dt_tmp <- merge(gene_list, dt_tmp)
			dt_tmp[,gene_symbol:=NULL]
			setkey(dt_tmp, 'geneset_name')
			dt_list[[i]] <- transpose(dt_tmp[, lapply(.SD, sum), by=geneset_name], make.names="geneset_name", keep.names="sample")
			start <- start + n_per_loop + 1
		}

		dt_test <- rbindlist(dt_list)
		fwrite(dt_test, file=paste0(file_out_prefix, "_", consequence, '_gene_lists.tsv'), sep='\t')
	}
}

create_exome_burden_file <- function(dt, consequence_categories, file_out_prefix="BPSCZ", burden_colname="coding_burden")
{
	dt_cons <- dt %>% filter(consequence_category %in% consequence_categories)
	dt_cons <- data.table(dt_cons)
	dt_cons[, c("consequence_category", "gene_symbol"):=NULL]
	dt_burden <- dt_cons[, lapply(.SD, sum)]
	dt_burden <- transpose(dt_burden, keep.names="sample")
	colnames(dt_burden)[2] <- burden_colname
	fwrite(dt_burden, file=paste0(file_out_prefix, "_burden.tsv"), sep='\t')
}
