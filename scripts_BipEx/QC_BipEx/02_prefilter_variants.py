import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Names of .mt files.
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
# Read in the hard calls matrix table.
mt = hl.read_matrix_table(MT_HARDCALLS)

# Read in the target intervals
TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list'
# Read in the padded target intervals (50bp padding)
PADDED_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets.interval_list'

# Low complexity regions in the data.
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/LCR-hs38.bed'

INITIAL_VARIANT_QC_FILE  = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_metrics.tsv'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

# Import the interval lists for the target intervals.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome='GRCh38')
# Import the interval lists for the padded target intervals.
padded_target_intervals = hl.import_locus_intervals(PADDED_TARGET_INTERVALS, reference_genome='GRCh38')
# Import the interval lists for the LCRs.
LCR_intervals = hl.import_locus_intervals(LCRs, reference_genome='GRCh38')

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals = ~hl.is_defined(padded_target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
in_LCR = mt.filter_rows(mt.in_LCR).count_rows()
not_in_target_intervals = mt.filter_rows(mt.not_in_target_intervals).count_rows()
not_in_padded_target_intervals = mt.filter_rows(mt.not_in_padded_target_intervals).count_rows()

print('n variants failing VQSR:')
pprint(fail_VQSR)
print('n variants in low complexity regions:')
pprint(in_LCR)
print('n variants not in target intervals:')
pprint(not_in_target_intervals)
print('n variants not in padded target intervals:')
pprint(not_in_padded_target_intervals)

# Export variant annotations, and variant keytable.
mt_rows = mt.rows()
mt_rows.select(mt_rows.fail_VQSR, mt_rows.in_LCR, mt_rows.not_in_padded_target_intervals).export(INITIAL_VARIANT_QC_FILE)
mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR | mt.not_in_padded_target_intervals, keep=False)

intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
mt = hl.filter_intervals(mt, intervals)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

mt_rows_filter = mt.rows().select().export(INITIAL_VARIANT_LIST)

n_variants = hl.import_table(INITIAL_VARIANT_LIST).count()

print('n variants after initial filter:')
print(n_variants)
