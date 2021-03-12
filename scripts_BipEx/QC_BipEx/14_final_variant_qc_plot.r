library(ggplot2)
library(ggsci)
library(data.table)

source("r_functions_and_parameters/r_options_BipEx.r")
source('r_functions_and_parameters/pretty_plotting.r')

save_figures <- TRUE

VARIANT_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.variants.tsv.bgz | gzcat'

df <- fread(VARIANT_QC_FILE, data.table=FALSE, header=TRUE, sep='\t')

df <- df %>% filter(qc.AF > 0 & qc.AF < 1)
df <- df %>% mutate(call_rate_diff = (case_qc.qc.call_rate - control_qc.qc.call_rate))

# call rate across all samples
create_pretty_hist(df, aes(x=qc.call_rate), threshold=T_variant_call_rate, x_label='Call Rate', xlim=c(0.9,1), save_figure=save_figures,
	file=paste0(PLOTS, '14_variant_call_rate'))
create_pretty_hist(df, aes(x=control_qc.qc.call_rate), threshold=T_variant_call_rate, x_label='Call Rate', xlim=c(0.9,1), save_figure=FALSE)
create_pretty_hist(df, aes(x=case_qc.qc.call_rate), threshold=T_variant_call_rate, x_label='Call Rate', xlim=c(0.9,1), save_figure=FALSE)

# cumulative call rate
p <- create_pretty_cumulative(df, aes(x=qc.call_rate, color='All samples'), x_label="Call Rate", threshold=T_variant_call_rate, 
    key_label='', xlim=c(0.9,1), save_figure=FALSE)
p <- p + stat_ecdf(aes(control_qc.qc.call_rate, color = 'Controls'), geom='line', pad=FALSE) + 
    stat_ecdf(aes(case_qc.qc.call_rate, color = 'Bipolar Disorder'), geom='line', pad=FALSE)

if(save_figures) {
	ggsave(paste0(PLOTS,'14_call_rate_cdf.jpg'), p, width=160, height=90, units='mm')
	ggsave(paste0(PLOTS,'14_call_rate_cdf.pdf'), p, width=160, height=90, units='mm')
}

# call rate difference (case rate - control rate)
create_pretty_hist(df, aes(x=call_rate_diff), x_label='Call Rate Difference',
    threshold=-T_absdiff, threshold_max=T_absdiff,
    xlim=c(-(1-T_variant_call_rate),(1-T_variant_call_rate)), save_figure=save_figures, file=paste0(PLOTS, '14_call_rate_diff'))

create_pretty_cumulative(df, aes(x=qc.p_value_hwe), x_label='p(HWE)',
    threshold=T_pHWE, xlim=c(0,1), save_figure=save_figures, file=paste0(PLOTS, '14_pHWE'))
