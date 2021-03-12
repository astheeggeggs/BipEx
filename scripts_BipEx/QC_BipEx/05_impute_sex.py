import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'

IMPUTESEX_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.ht'
IMPUTESEX_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.tsv'
Y_NCALLED = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_ycalled.tsv'

INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'
PRUNED_CHRX_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/04_chrX.prune.in'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

ht_pruned_chrx_variants = ht_pruned_chrx_variants.annotate(**hl.parse_variant(ht_pruned_chrx_variants.f0, reference_genome='GRCh38'))
ht_pruned_chrx_variants = ht_pruned_chrx_variants.key_by(ht_pruned_chrx_variants.locus, ht_pruned_chrx_variants.alleles)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
# Want to change this to reflect the dataset that I have.
mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# Determine non-missing allele count on the y.
mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)
