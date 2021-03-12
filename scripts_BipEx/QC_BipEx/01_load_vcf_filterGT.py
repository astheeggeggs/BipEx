import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# The following was run using these jars and zips (to enable the new version of split_multi).
# gs://hail-common/builds/0.2/jars/hail-0.2-7a280b932abc959f7accf11124f3255038f48c95-Spark-2.4.0.jar
# gs://hail-common/builds/0.2/python/hail-0.2-7a280b932abc959f7accf11124f3255038f48c95.zip

RAW_MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes.mt'
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'

# Remove multiallelics with 100 or more alleles...
mt = hl.read_matrix_table(RAW_MT)

# Count before splitting multi-allelics.
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

# Read in the target intervals
TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list'
# Import the interval lists for the target intervals.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome='GRCh38')

mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
not_in_target_intervals = mt.filter_rows(mt.not_in_target_intervals).count_rows()

print('')
print('n variants not in target intervals:')
pprint(not_in_target_intervals)

# Read in the padded target intervals (50bp padding)
PADDED_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets.interval_list'
# Import the interval lists for the padded target intervals.
padded_target_intervals = hl.import_locus_intervals(PADDED_TARGET_INTERVALS, reference_genome='GRCh38')

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(not_in_padded_target_intervals = ~hl.is_defined(padded_target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
not_in_padded_target_intervals = mt.filter_rows(mt.not_in_padded_target_intervals).count_rows()

print('')
print('n variants not in padded target intervals:')
pprint(not_in_padded_target_intervals)

# Need to do the following for testing.
mt = mt.filter_rows(mt.alleles.length() <= 6)

n = mt.count_rows()

pprint('')
pprint('n variants not more than 6 alleles:')
print(n)

mt = hl.split_multi_hts(mt)

mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                # ((mt.AD[0] / mt.DP) < 0.8) | # Has to be removed because allele depth no longer defined for hom ref calls.
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)

mt.write(MT, overwrite=True)
mt = mt.checkpoint(MT, overwrite=True)
mt = hl.read_matrix_table(MT)
mt.select_entries(mt.GT).repartition(512).write(MT_HARDCALLS, overwrite=True)

mt = hl.read_matrix_table(MT_HARDCALLS)
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])
