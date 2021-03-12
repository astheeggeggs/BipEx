import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'
ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'
URV_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'
SEXCHECK_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_sexcheck.remove.sample_list'
IBD_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list'

ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST,
    types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(locus=ht_initial_variants.locus, alleles=ht_initial_variants.alleles)
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')

# Annotate variants with counts from non-psychiatric version of Gnomad.
# Fill in variants not in Gnomad variant list.
# Annotate variants with LoF/damaging missense annotation.

mt = hl.read_matrix_table(MT)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
# Drop some fields that are not needed.
mt = mt.drop('a_index', 'qual', 'rsid', 'info', 'filters', 'was_split')

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
constraint_annotations = hl.read_table(ANNOTATION_TABLE)

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_rows(constraint = constraint_annotations[mt.row_key])
mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
mt = mt.filter_rows((mt.is_singleton) & (~mt.constraint.inGnomAD_nonpsych))

mt = mt.annotate_cols(n_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "ptv")),
                      n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "damaging_missense")),
                      n_URV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "other_missense")),
                      n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "synonymous")),
                      n_URV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "non_coding")))

mt.cols().flatten().export(URV_FILE)
