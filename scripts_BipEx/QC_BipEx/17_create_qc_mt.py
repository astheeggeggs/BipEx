import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'

PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
IMPUTESEX_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.ht'

ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

FINAL_SAMPLE_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.keep.sample_list'
FINAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.keep.variant_list'

FINAL_PRUNED_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/16_prune.final_qc.keep.variant_list'

QC_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.mt'
QC_HARDCALLS_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls.mt'

FINAL_VARIANT_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/17_final_qc.variants.tsv.bgz'
FINAL_SAMPLE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0')
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

ht_final_pruned_variants = hl.import_table(FINAL_PRUNED_VARIANTS, no_header=True)
ht_final_pruned_variants = ht_final_pruned_variants.annotate(**hl.parse_variant(ht_final_pruned_variants.f0, reference_genome='GRCh38'))
ht_final_pruned_variants = ht_final_pruned_variants.key_by(ht_final_pruned_variants.locus, ht_final_pruned_variants.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)
annotation_annotations = hl.read_table(ANNOTATION_TABLE)

mt = hl.read_matrix_table(MT)
mt = mt.drop('a_index', 'qual', 'info', 'filters', 'was_split')

mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

mt = mt.annotate_cols(phenotype = sample_annotations[mt.col_key])
mt = mt.annotate_cols(imputesex = impute_sex_annotations[mt.col_key])
mt = mt.annotate_rows(annotation = annotation_annotations[mt.row_key])

mt = hl.variant_qc(mt, name = 'qc')

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.het_freq_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).het_freq_hwe)))
)

mt = mt.annotate_rows(
	annotation = mt.annotation.annotate(
		info = mt.annotation.info.annotate(
			AC = mt.annotation.info.AC[mt.annotation.a_index-1],
			AF = mt.annotation.info.AF[mt.annotation.a_index-1],)
		)
	)

mt = hl.sample_qc(mt)

mt_pca = mt.filter_rows(hl.is_defined(ht_final_pruned_variants[mt.row_key]))

pca_output = hl.hwe_normalized_pca(mt_pca.GT, k=10)
pca_output = pca_output[1].key_by('s')
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

mt = mt.annotate_cols(pca = pca_output[mt.s])

n = mt.count()

print('')
print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt.write(QC_MT, overwrite=True)
mt.select_entries(mt.GT).repartition(512).write(QC_HARDCALLS_MT, overwrite=True)
