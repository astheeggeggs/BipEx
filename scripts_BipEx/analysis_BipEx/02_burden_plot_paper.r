library(data.table)
source("../QC_BipEx/r_functions_and_parameters/pretty_plotting.r")

dir.create("../../analysis_plots/forest_plots/point_range_plots", recursive=TRUE, showWarnings=FALSE)
load("forest_plot_Rdata_files/BP_including_BPSCZ.Rdata")

# Make into a data.table for plotting
init <- TRUE
for(pheno in names(forest_plots)) {
	for (test in c("lin", "logit")) {
		if(init) {
			forest_plots_coding_burden_dt <- data.table(forest_plots_coding_burden[[pheno]][[test]], phenotype=pheno, test=test)
			init <- FALSE
		} else {
			forest_plots_coding_burden_dt <- rbind(forest_plots_coding_burden_dt, data.table(forest_plots_coding_burden[[pheno]][[test]], phenotype=pheno, test=test))
		}
	}
}

add_descriptions <- function(dt)
{
	return(
		dt %>% mutate(
			phenotype_description = case_when(
				phenotype=='is_BP' ~ "Bipolar Disorder",
				phenotype=='is_BP_including_BPSCZ' ~ "Bipolar Disorder (including Schizoaffective)",
				phenotype=='is_BP1' ~ "Bipolar Disorder 1",
				phenotype=='is_BP2' ~ "Bipolar Disorder 2",
				phenotype=='is_BPNOS' ~ "Bipolar Disorder NOS",
				phenotype=='is_BPSCZ' ~ "Schizoaffective",
				phenotype=='is_BPPSY' ~ "Bipolar Disorder, with Psychosis",
				phenotype=='is_BP_no_PSY' ~ "Bipolar Disorder, without Psychosis"
			),

			pLI = case_when(
				grepl("pli_09\\.", label) ~ "Not in gnomAD\npLI > 0.9",
				grepl("pli_0995\\.", label) ~ "Not in gnomAD\npLI > 0.995",
				grepl("gnom_non_psych\\.", label) ~ "Not in gnomAD",
				TRUE ~ "All"),

			label_description = case_when(
				grepl("URV_damaging_missense", label) ~ "Damaging missense",
				grepl("PTV", label) ~ "PTV",
				grepl("other_missense", label) ~ "Other missense",
				grepl("MPC_2_missense", label) ~ "Missense MPC > 2",
				grepl("synonymous", label) ~ "Synonymous",
				grepl("non_coding", label) ~ "Non-coding",
				TRUE ~ ""
			)
		)
	)
}

# All genes, constrained genes (pLI > 0.9), constrained genes (pLI > 0.995).
# Include non-coding, synonymous, other missense, damaging missense, MPC > 2.

forest_plots_coding_burden_dt <- add_descriptions(forest_plots_coding_burden_dt)

# Select the consequence classes to include
# Create the names for them in the legend.

# Include the relevant y-line.

# Loop and generate.

dt <- forest_plots_coding_burden_dt
dt <- dt %>% filter(
	grepl('PTV', label) |
	grepl('MPC_2_missense', label) |
	grepl('URV_damaging_missense', label) |
	grepl('other_missense', label) |
	grepl('synonymous', label)
) %>% filter(pLI != "Not in gnomAD\npLI > 0.995") %>% 
	mutate(pLI = factor(pLI,
		levels=c("Not in gnomAD\npLI > 0.9",
			"Not in gnomAD",
			"All")))


i <- 1

pheno_id <- c("is_BP", "is_BP1", "is_BP2")
pheno_description <- c("Bipolar Disorder", "Bipolar Disorder 1", "Bipolar Disorder 2")
pdf("pointrange_plots_gnomAD.pdf", width=5, height=3.5)

for(pheno in pheno_description) {
	for(current_test in c("logit", "lin")) {
		if(current_test == "logit") {
			y_label <- "log(odds ratio)"
			ylim <- c(-0.05,0.1)
			threshold <- 0
			dt_current <- dt %>% filter(test == current_test) %>% mutate(mean = log10(mean), lower=log10(lower), upper=log10(upper))
		} else {
			y_label <- "Excess variants"
			ylim <- c(-0.25,0.25)
			threshold <- 0
			dt_current <- dt %>% filter(test == current_test)
		}
		create_pretty_pointrange(
			dt_current %>% filter(phenotype_description == pheno),
			"mean", "lower", "upper", title=pheno, break_on="pLI",
			colors="label_description", title.hjust=0, y_label=y_label,
			print_p = TRUE, title_size=18, threshold=threshold, save_figure=TRUE,
			colour_levels = c("PTV", "Missense MPC > 2", "Damaging missense", "Other missense", "Synonymous"),
			ylim=ylim,
			manual_colours = rev(c("#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D"))
		)
	}
	i <- i+1
}

dev.off()


i <- 1
pLI_id <- unique(dt$pLI)
pLI_title <- gsub("\\\n", ", ", pLI_id)
file_names <- gsub("[ >\\.\\\n]+", "_", pLI_id)

dt <- dt %>% filter(phenotype_description %in% c("Bipolar Disorder", "Bipolar Disorder 1", "Bipolar Disorder 2"))
dt$label_description[which(dt$label_description == "Missense MPC > 2")] <- "Missense\nMPC > 2"
dt$label_description[which(dt$label_description == "Damaging missense")] <- "Damaging\nmissense"

pdf("pointrange_plots_classes.pdf", width=7.5, height=3.5)

for(conservation in unique(dt$pLI)) {
	for(current_test in c("logit", "lin")) {

		if(current_test == "logit") {
			y_label <- "log(odds ratio)"
			threshold <- 0
			dt_current <- dt %>% filter(test == current_test) %>% mutate(mean = log10(mean), lower=log10(lower), upper=log10(upper))
		} else {
			y_label <- "Excess variants"
			threshold <- 0
			dt_current <- dt %>% filter(test == current_test)
		}

		create_pretty_pointrange(
			dt_current %>% filter(pLI == conservation),
			"mean", "lower", "upper", title=pLI_title[i], break_on="label_description",
			colors="phenotype_description", title.hjust=0, y_label=y_label,
			print_p = TRUE, title_size=18, threshold=threshold, save_figure=TRUE,
			colour_levels = c("Bipolar Disorder", "Bipolar Disorder 1", "Bipolar Disorder 2"),
			break_on_levels = c("PTV", "Missense\nMPC > 2", "Damaging\nmissense", "Other missense", "Synonymous")
		)
	}
	i <- i+1
}
dev.off()

