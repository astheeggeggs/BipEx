import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
IBD_OUTPUT = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.tsv'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
PRUNED_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/04_prune.keep.variant_list'

INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh38'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s]).repartition(128).persist()

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

ibd_table = hl.identity_by_descent(mt, min=0.1)
ibd_table = ibd_table.flatten()
ibd_table.export(IBD_OUTPUT)
