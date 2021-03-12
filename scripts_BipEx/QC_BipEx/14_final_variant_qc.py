import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
IMPUTESEX_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.ht'

# NOTE: IF THINGS LOOKED CLEAN FROM THE PREVIOUS PCA PLOT,
# DEFINE PCA_LIST TO BE THE SAME AS THE STRICT EUROPEAN SAMPLE LIST.
# OTHERWISE, RESTRICT.

# Also, need to ensure that the AJ samples used to train are also removed.

PCA_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european.sample_list'
# These are the AJs that are classified using the AJs from the US samples.
AJ_CLASSIFY_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_aj_classify.sample_list'
# These are the samples that we identify as AJ in the US samples.
AJ_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_aj.sample_list'

IBD_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list'

VARIANT_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.variants.tsv.bgz'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)

# This is empty - uncomment if it's not.
# ht_aj_samples = hl.import_table(AJ_CLASSIFY_LIST, no_header=True, key='f0')
ht_aj_samples = hl.import_table(AJ_LIST, no_header=True, key='f0')
ht_pca_samples = hl.import_table(PCA_LIST, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')

ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})

ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

mt = hl.read_matrix_table(MT)

mt = mt.filter_cols(hl.is_defined(ht_pca_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_aj_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))

mt = mt.annotate_cols(phenotype = sample_annotations[mt.col_key])
mt = mt.annotate_cols(imputesex = impute_sex_annotations[mt.col_key])

mt = hl.variant_qc(mt, name='qc')

mt = mt.annotate_rows(
	qc=mt.qc.annotate(AC=mt.qc.AC[1],
	AF=mt.qc.AF[1],
	homozygote_count=mt.qc.homozygote_count[1]))

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.p_value_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).p_value)))
)

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.het_freq_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).het_freq_hwe)))
)

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

# Want to compare the average quality of variants in cases and controls, and compare them.
# We also want to compare the quality of the variants across the different locations.

mt_control = mt.filter_cols(mt.phenotype.PHENOTYPE_COARSE == "Control")
mt_control = hl.variant_qc(mt_control, name = 'qc')

mt_control = mt_control.annotate_rows(
	qc=mt_control.qc.annotate(AC=mt_control.qc.AC[1],
	AF=mt_control.qc.AF[1],
	homozygote_count=mt_control.qc.homozygote_count[1]))

mt_control = mt_control.annotate_rows(qc = mt_control.qc.annotate(p_value_hwe = hl.case()
	.when(mt_control.locus.in_autosome(), mt_control.qc.het_freq_hwe)
	.default(hl.agg.filter(mt_control.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt_control.GT).het_freq_hwe)))
)

n = mt_control.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_case = mt.filter_cols(mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder")
mt_case = hl.variant_qc(mt_case, name='qc')

mt_case = mt_case.annotate_rows(
	qc=mt_case.qc.annotate(AC=mt_case.qc.AC[1],
	AF=mt_case.qc.AF[1],
	homozygote_count=mt_case.qc.homozygote_count[1]))

mt_case = mt_case.annotate_rows(qc = mt_case.qc.annotate(p_value_hwe = hl.case()
	.when(mt_case.locus.in_autosome(), mt_case.qc.het_freq_hwe)
	.default(hl.agg.filter(mt_case.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt_case.GT).het_freq_hwe)))
)

n = mt_case.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

# Finally, annotate the larger matrix table with this information, and export.
mt_control_rows = mt_control.rows().select('qc')
mt = mt.annotate_rows(control_qc = mt_control_rows[mt.row_key])

mt_case_rows = mt_case.rows().select('qc')
mt = mt.annotate_rows(case_qc = mt_case_rows[mt.row_key])

mt.rows().select('qc', 'control_qc', 'case_qc').flatten().export(VARIANT_QC_FILE)

