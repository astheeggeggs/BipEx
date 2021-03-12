import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.mt' 
SAMPLE_QC_IN_TARGET = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/18_final_qc.after_in_target.samples.tsv'

mt = hl.read_matrix_table(QC_MT)

mt_het = mt.filter_entries(mt.GT.is_het())
het_struct = mt_het.aggregate_entries(hl.struct(global_ADhet=hl.agg.stats(mt_het.AD[1]/mt_het.DP),
	global_ADhet_20 = hl.agg.mean(mt_het.AD[1]/mt_het.DP < 0.2),
	global_ADhet_25 = hl.agg.mean(mt_het.AD[1]/mt_het.DP < 0.25),
	global_ADhet_30 = hl.agg.mean(mt_het.AD[1]/mt_het.DP < 0.30),
	global_ADhet_35 = hl.agg.mean(mt_het.AD[1]/mt_het.DP < 0.35),
	global_ADhom=hl.agg.stats((mt_het.AD[1] + mt_het.AD[1])/mt_het.DP)))

print(het_struct)

mt_hom_var = mt.filter_entries(mt.GT.is_hom_var())
hom_struct = mt_hom_var.aggregate_entries(hl.struct(global_ADhet=hl.agg.stats(mt_hom_var.AD[1]/mt_hom_var.DP),
	global_ADhom=hl.agg.stats((mt_hom_var.AD[0] + mt_hom_var.AD[1])/mt_hom_var.DP),
	global_ADhom_80 = hl.agg.mean((mt_hom_var.AD[1] + mt_hom_var.AD[1])/mt_hom_var.DP < 0.8),
	global_ADhom_85 = hl.agg.mean((mt_hom_var.AD[1] + mt_hom_var.AD[1])/mt_hom_var.DP < 0.85),
	global_ADhom_90 = hl.agg.mean((mt_hom_var.AD[1] + mt_hom_var.AD[1])/mt_hom_var.DP < 0.9),
	global_ADhom_95 = hl.agg.mean((mt_hom_var.AD[1] + mt_hom_var.AD[1])/mt_hom_var.DP < 0.95)))
print(hom_struct)

# Also, wish to examine the Ti/Tv ratio within the calling intervals (excluding the padding).

mt = mt.filter_rows(~mt.not_in_target_intervals)
mt = hl.sample_qc(mt, name='sample_qc_in_target')
mt.cols().select("imputesex", "sample_qc", "sample_qc_in_target").flatten().export(SAMPLE_QC_IN_TARGET)
