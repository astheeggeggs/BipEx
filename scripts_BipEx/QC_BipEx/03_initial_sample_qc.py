import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Create sample QC metrics restricted and not restricted to the ICE intervals.

# Inputs
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

# Outputs
INITIAL_SAMPLE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv'

variants_to_filter = hl.import_table(INITIAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)

mt = hl.read_matrix_table(MT)
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])

n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

mt = hl.sample_qc(mt, name='qc_padded_ice')

TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list'
# Import the interval lists for the LCRs.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome='GRCh38')
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)

n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

mt = hl.sample_qc(mt, name='qc_ice')
mt.cols().select('phenotype', 'qc_padded_ice', 'qc_ice').flatten().export(output=INITIAL_SAMPLE_QC_FILE)

