import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Names of .mt files.
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.hardcalls.mt'
# Read in the hard calls matrix table.
mt = hl.read_matrix_table(MT_HARDCALLS)

# Read in the target intervals
TARGET_INTERVALS = 'gs://dalio_bipolar_w1_w2_hail_02/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.fixed.interval_list'

# Low complexity regions in the data.
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/inputs_low_complexity_regions_b37.interval_list'

INITIAL_VARIANT_QC_FILE  = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_metrics_b37_callset.tsv'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_b37_callset.keep.variant_list'
INITIAL_VARIANT_AUTO_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_b37_callset.keep.autosome.variant_list'

# Import the interval lists for the LCRs.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS)
LCR_intervals = hl.import_locus_intervals(LCRs)

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
in_LCR = mt.filter_rows(mt.in_LCR).count_rows()
not_in_target_intervals = mt.filter_rows(mt.not_in_target_intervals).count_rows()

print('n variants failing VQSR:')
pprint(fail_VQSR)
print('n variants in low complexity regions:')
pprint(in_LCR)
print('n variants not in target intervals:')
pprint(not_in_target_intervals)

# Export variant annotations, and variant keytable.
mt_rows = mt.rows()
mt_rows.select(mt_rows.fail_VQSR, mt_rows.in_LCR, mt_rows.not_in_target_intervals).export(INITIAL_VARIANT_QC_FILE)

mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR | mt.not_in_target_intervals, keep=False)

# Filter to variants in the autosomes.
mt_auto = hl.filter_intervals(mt,
	[hl.parse_locus_interval('1:START-22:END')])

intervals = [hl.parse_locus_interval(x) for x in ['1:START-22:END', 'X:START-X:END', 'Y:START-Y:END']]
mt = hl.filter_intervals(mt, intervals)

mt_auto_rows_filter = mt_auto.rows().select().export(INITIAL_VARIANT_AUTO_LIST)
mt_rows_filter = mt.rows().select().export(INITIAL_VARIANT_LIST)

n_variants = hl.import_table(INITIAL_VARIANT_LIST).count()

print('n variants after initial filter:')
print(n_variants)

n_variants = hl.import_table(INITIAL_VARIANT_AUTO_LIST).count()

print('n variants after initial filter and restricting to autosomes:')
print(n_variants)
