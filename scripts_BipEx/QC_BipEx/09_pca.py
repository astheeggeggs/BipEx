import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
PCA_SCORES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/09_pca_scores.tsv'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'
PRUNED_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/04_prune.keep.variant_list'
IBD_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list'

mt = hl.read_matrix_table(MT_HARDCALLS)
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh38'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s]).repartition(128).persist()

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

pca_output = hl.hwe_normalized_pca(mt.GT, k=10)
pca_output = pca_output[1].key_by('s')

pca_output = pca_output.annotate(phenotype = sample_annotations[pca_output.s]).repartition(128).persist()

pca_output = pca_output.annotate(PC1 = pca_output.scores[0],
	PC2 = pca_output.scores[1], PC3 = pca_output.scores[2],
	PC4 = pca_output.scores[3], PC5 = pca_output.scores[4],
	PC6 = pca_output.scores[5], PC7 = pca_output.scores[6],
	PC8 = pca_output.scores[7], PC9 = pca_output.scores[8],
	PC10 = pca_output.scores[9])

pca_output.flatten().export(output=PCA_SCORES)
