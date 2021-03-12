import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
MT_1KG = 'gs://hail-datasets-hail-data/1000_Genomes_autosomes.phase_3.GRCh38.mt'
PCA_SCORES_EUR = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.european.1kg.tsv'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
POPULATIONS_1KG = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/samples_1kg.ht'
INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'
PRUNED_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/04_prune.keep.variant_list'

IBD_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list'
EUROPEAN_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european.sample_list'

mt_1kg = hl.read_matrix_table(MT_1KG)
mt_1kg = hl.split_multi_hts(mt_1kg)
# This is to enable a join later.
mt_1kg = mt_1kg.select_entries("GT")

populations_1kg = hl.read_table(POPULATIONS_1KG)
mt_1kg = mt_1kg.annotate_cols(pop_1kg = populations_1kg[mt_1kg.s])
mt_1kg = mt_1kg.filter_cols(mt_1kg.pop_1kg.super_population == "EUR")
mt_1kg = mt_1kg.select_cols()

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh38'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_european_samples = hl.import_table(EUROPEAN_SAMPLES, no_header=True, key='f0')

mt = hl.read_matrix_table(MT_HARDCALLS)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_cols(hl.is_defined(ht_european_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))

mt = mt.union_cols(mt_1kg)

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

pca_output = hl.hwe_normalized_pca(mt.GT, k=10)
pca_output = pca_output[1].key_by('s')

sample_annotations = hl.read_table(PHENOTYPES_TABLE)

pca_output = pca_output.annotate(phenotype = sample_annotations[pca_output.s])
pca_output = pca_output.annotate(
	PC1 = pca_output.scores[0],
	PC2 = pca_output.scores[1],
	PC3 = pca_output.scores[2],
	PC4 = pca_output.scores[3],
	PC5 = pca_output.scores[4],
	PC6 = pca_output.scores[5],
	PC7 = pca_output.scores[6],
	PC8 = pca_output.scores[7],
	PC9 = pca_output.scores[8],
	PC10 = pca_output.scores[9])

pca_output = pca_output.annotate(pop_1kg = populations_1kg[pca_output.s])

pca_output = pca_output.annotate(PHENOTYPE_COARSE = hl.case().when(hl.is_defined(pca_output.phenotype.PHENOTYPE_COARSE), pca_output.phenotype.PHENOTYPE_COARSE)
	.default("1KG"),
	PI = hl.case().when(hl.is_defined(pca_output.phenotype.PI), pca_output.phenotype.PI)
	.default("1KG"),
	LOCATION = hl.case().when(hl.is_defined(pca_output.phenotype.LOCATION), pca_output.phenotype.LOCATION)
	.default("1KG"),
	PHENOTYPE_FINE = hl.case().when(hl.is_defined(pca_output.phenotype.PHENOTYPE_FINE), pca_output.phenotype.PHENOTYPE_FINE)
	.default("1KG"),
	PROJECT_OR_COHORT = hl.case().when(hl.is_defined(pca_output.phenotype.PROJECT_OR_COHORT), pca_output.phenotype.PROJECT_OR_COHORT)
	.default("1KG"),
	SUPER_POPULATION = hl.case().when(hl.is_defined(pca_output.pop_1kg.super_population), pca_output.pop_1kg.super_population)
	.default(pca_output.phenotype.PHENOTYPE_COARSE),
	POPULATION = hl.case().when(hl.is_defined(pca_output.pop_1kg.population), pca_output.pop_1kg.population)
	.default(pca_output.phenotype.PHENOTYPE_COARSE)
).repartition(128).persist()

pca_output.flatten().export(output=PCA_SCORES_EUR)
