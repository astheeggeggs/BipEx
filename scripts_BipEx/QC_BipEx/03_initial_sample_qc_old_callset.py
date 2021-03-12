import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Inputs
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.mt'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_b37_callset.keep.variant_list'
INITIAL_VARIANT_AUTO_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_b37_callset.keep.autosome.variant_list'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes_old_callset.ht'

# Outputs
INITIAL_SAMPLE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc_b37_callset.tsv'
INITIAL_SAMPLE_QC_FILE_INV_REMOVED = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc_b37_callset_inv_removed.tsv'

variants_to_filter = hl.import_table(INITIAL_VARIANT_AUTO_LIST,
	types={'locus':hl.tlocus(), 'alleles':hl.tarray(hl.tstr)})
variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)

pprint(sample_annotations.describe())

mt = hl.read_matrix_table(MT)
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt_invariant_included = hl.sample_qc(mt, name='qc')

mt_invariant_included.cols().select('phenotype', 'qc').flatten().export(output=INITIAL_SAMPLE_QC_FILE)

mt = hl.variant_qc(mt, name='qc')
mt_invariant_removed = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))
mt_invariant_removed = hl.sample_qc(mt_invariant_removed, name='qc_sample')

mt_invariant_removed.cols().select('phenotype', 'qc_sample').flatten().export(output=INITIAL_SAMPLE_QC_FILE_INV_REMOVED)
