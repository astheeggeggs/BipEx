import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Annotate variants with counts from non-psychiatric version of Gnomad.
# Fill in variants not in Gnomad variant list.
# Annotate variants with LoF/damaging missense annotation.

MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes.mt'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'
ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

URV_NOT_IN_GNOMAD_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/20_URVs_not_in_gnomAD.tsv'
NOT_IN_GNOMAD_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/20_not_in_gnomAD.tsv'

mt = hl.read_matrix_table(MT)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
constraint_annotations = hl.read_table(ANNOTATION_TABLE)

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_rows(constraint = constraint_annotations[mt.row_key])
mt = mt.filter_rows(~mt.constraint.inGnomAD_nonpsych)

mt = mt.annotate_cols(n_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_coding_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "ptv")),
                      n_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "damaging_missense")),
                      n_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "other_missense")),
                      n_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "synonymous")),
                      n_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "non_coding")))

mt.cols().flatten().export(NOT_IN_GNOMAD_FILE)

mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
mt = mt.filter_rows(mt.is_singleton)

mt = mt.annotate_cols(n_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "ptv")),
                      n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "damaging_missense")),
                      n_URV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "other_missense")),
                      n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "synonymous")),
                      n_URV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "non_coding")))

mt.cols().flatten().export(URV_NOT_IN_GNOMAD_FILE)
