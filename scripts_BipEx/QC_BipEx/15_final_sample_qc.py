import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
IMPUTESEX_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.ht'

SEXCHECK_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_sexcheck.remove.sample_list'

# NOTE: IF THINGS LOOKED CLEAN FROM THE PREVIOUS PCA PLOT,
# DEFINE PCA_LIST TO BE THE SAME AS THE STRICT EUROPEAN SAMPLE LIST.
# OTHERWISE, RESTRICT.

# HERE, IT WAS CLEAN, SO I DEFINE IT AS THE LIST OF STRICTLY DEFINED EUROPEANS.
PCA_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european.sample_list'
AJ_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_aj.sample_list'
# AJ_CLASSIFY_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_AJ_classify.sample_list'
IBD_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list'

INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'
FINAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.keep.variant_list'

SAMPLE_BEFORE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.after.samples.tsv'

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)

# This is empty - uncomment if it's not.
# ht_aj_samples = hl.import_table(AJ_CLASSIFY_LIST, no_header=True, key='f0')
ht_aj_samples = hl.import_table(AJ_LIST, no_header=True, key='f0')
ht_pca_samples = hl.import_table(PCA_LIST, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_sex_check_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})

ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

ht_final_variants = hl.import_table(FINAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})

ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

mt_before = hl.read_matrix_table(MT)
mt_before = mt_before.filter_cols(hl.is_defined(ht_pca_samples[mt_before.col_key]))

mt_before = mt_before.filter_cols(~hl.is_defined(ht_aj_samples[mt_before.col_key]))
mt_before = mt_before.filter_cols(~hl.is_defined(ht_ibd_samples[mt_before.col_key]))
mt_before = mt_before.filter_cols(~hl.is_defined(ht_sex_check_samples[mt_before.col_key]))
mt_before = mt_before.filter_rows(hl.is_defined(ht_initial_variants[mt_before.row_key]))

mt_before = mt_before.annotate_cols(phenotype = sample_annotations[mt_before.col_key])
mt_before = mt_before.annotate_cols(imputesex = impute_sex_annotations[mt_before.col_key])

mt_before = hl.variant_qc(mt_before, name = 'qc')

mt_before = mt_before.annotate_rows(
	qc = mt_before.qc.annotate(AC=mt_before.qc.AC[1],
	AF = mt_before.qc.AF[1],
	homozygote_count = mt_before.qc.homozygote_count[1]))

mt_before = mt_before.filter_rows((mt_before.qc.AF > 0) & (mt_before.qc.AF < 1))
mt_before = hl.sample_qc(mt_before)

n = mt_before.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_before = mt_before.annotate_cols(sex = hl.case()
	.when(mt_before.imputesex.impute_sex.is_female, "Female")
	.default("Male"))

mt_after = mt_before.filter_rows(hl.is_defined(ht_final_variants[mt_before.row_key]))
mt_after = hl.sample_qc(mt_after)

n = mt_after.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_before.cols().select("sex", "phenotype", "sample_qc").flatten().export(SAMPLE_BEFORE_QC_FILE)
mt_after.cols().select("sex", "phenotype", "sample_qc").flatten().export(SAMPLE_AFTER_QC_FILE)
