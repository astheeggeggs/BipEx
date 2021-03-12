import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Names of .mt files.
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'
# Read in the hard calls matrix table.
mt = hl.read_matrix_table(MT)

# Read in the target intervals
PADDING_0_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list'
# Read in the padded target intervals (50bp padding)
PADDING_50_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets.interval_list'
# Read in the padded target intervals (100bp padding)
PADDING_100_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets_100.interval_list'
# Read in the padded target intervals (150bp padding)
PADDING_150_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets_150.interval_list'

# Low complexity regions in the data.
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/LCR-hs38.bed'

PADDED_150_INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/padding.150_padding.variant_list'
PADDED_100_INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/padding.100_padding.variant_list'
PADDED_50_INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/padding.50_padding.variant_list'
PADDED_0_INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/padding.0_padding.variant_list'

PADDING_METRICS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/padding_sample_metrics_ICE.tsv.gz'

# Import the interval lists for the target intervals.
target_intervals = hl.import_locus_intervals(PADDING_0_TARGET_INTERVALS, reference_genome='GRCh38')

# Import the interval lists for the padded target intervals.
padded_target_intervals_50 = hl.import_locus_intervals(PADDING_50_TARGET_INTERVALS, reference_genome='GRCh38')
padded_target_intervals_100 = hl.import_locus_intervals(PADDING_100_TARGET_INTERVALS, reference_genome='GRCh38')
padded_target_intervals_150 = hl.import_locus_intervals(PADDING_150_TARGET_INTERVALS, reference_genome='GRCh38')

# Import the interval lists for the LCRs.
LCR_intervals = hl.import_locus_intervals(LCRs, reference_genome='GRCh38')

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals_50 = ~hl.is_defined(padded_target_intervals_50[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals_100 = ~hl.is_defined(padded_target_intervals_100[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals_150 = ~hl.is_defined(padded_target_intervals_150[mt.locus]))

# Export variant annotations, and variant keytable.
mt_rows = mt.rows()
mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR | mt.not_in_padded_target_intervals_150, keep=False)

intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
mt = hl.filter_intervals(mt, intervals)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

mt_rows_filter = mt.rows().select().export(PADDED_150_INITIAL_VARIANT_LIST)

n_variants = hl.import_table(PADDED_150_INITIAL_VARIANT_LIST).count()
print('n variants after initial filter:')
print(n_variants)

mt = hl.sample_qc(mt, name='qc_150')

mt = mt.filter_rows(mt.not_in_padded_target_intervals_100, keep=False)
mt_rows_filter = mt.rows().select().export(PADDED_100_INITIAL_VARIANT_LIST)
n_variants = hl.import_table(PADDED_100_INITIAL_VARIANT_LIST).count()
print('n variants after initial filter:')
print(n_variants)

mt = hl.sample_qc(mt, name='qc_100')

mt = mt.filter_rows(mt.not_in_padded_target_intervals_50, keep=False)
mt_rows_filter = mt.rows().select().export(PADDED_50_INITIAL_VARIANT_LIST)
n_variants = hl.import_table(PADDED_50_INITIAL_VARIANT_LIST).count()
print('n variants after initial filter:')
print(n_variants)

mt = hl.sample_qc(mt, name='qc_50')

mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)
mt_rows_filter = mt.rows().select().export(PADDED_0_INITIAL_VARIANT_LIST)
n_variants = hl.import_table(PADDED_0_INITIAL_VARIANT_LIST).count()
print('n variants after initial filter:')
print(n_variants)

mt = hl.sample_qc(mt, name='qc_0')

mt_rows_filter = mt.cols().flatten().export(PADDING_METRICS)