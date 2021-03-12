library(data.table)

df_old <- fread(cmd='gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc_b37_callset.tsv')
df_new <- fread(cmd='gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv')

df <- merge(df_old, df_new, by='s')

plot(df$qc.n_hom_ref.x, df$qc.n_hom_ref.y)
plot(df$qc.n_het.x, df$qc.n_het.y)
plot(df$qc.n_hom_var.x, df$qc.n_hom_var.y)

plot(df$qc.call_rate.x, df$qc.call_rate.y)
plot(df$qc.dp_stats.mean.x, df$qc.dp_stats.mean.y)
plot(df$qc.gq_stats.mean.x, df$qc.gq_stats.mean.y)

plot(df$qc.n_called.x, df$qc.n_called.y)
plot(df$qc.n_not_called.x, df$qc.n_not_called.y)
plot(df$qc.n_filtered.x, df$qc.n_filtered.y)

plot(df$qc.n_insertion.x, df$qc.n_insertion.y)
plot(df$qc.n_deletion.x, df$qc.n_deletion.y)

plot(df$qc.n_singleton.x, df$qc.n_singleton.y)
plot(df$qc.n_non_ref.x, df$qc.n_non_ref.y)
plot(df$qc.n_star.x, df$qc.n_star.y)
plot(df$qc.n_transition.x, df$qc.n_transition.y)
plot(df$qc.n_transversion.x, df$qc.n_transversion.y)

plot(df$qc.r_het_hom_var.x, df$qc.r_het_hom_var.y)
plot(df$qc.r_insertion_deletion.x, df$qc.r_insertion_deletion.y)
plot(df$qc.r_ti_tv.x, df$qc.r_ti_tv.y)