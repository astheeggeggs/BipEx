import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

ht = hl.read_table(ANNOTATION_TABLE)
ht = ht.annotate(type = hl.case().when((ht.alleles[0].length() == 1) | (ht.alleles[1].length() == 1), "SNP")
	.when(ht.alleles[0].length() < ht.alleles[1].length(), "Insertion")
	.when(ht.alleles[0].length() > ht.alleles[1].length(), "Deletion")
	.default("No type"))

n = ht.count()
n_inGnomAD_nonpsych = ht.aggregate(hl.agg.counter(ht.inGnomAD_nonpsych))
n_have_gene = ht.aggregate(hl.agg.counter(hl.is_defined(ht.vep.worst_csq_for_variant_canonical.gene_symbol)))
n_consequence_category = ht.aggregate(hl.agg.counter(ht.consequence_category))
n_most_severe_consequence = ht.aggregate(hl.agg.counter(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence))

print('n variants:')
print(n)

print('gnomAD_nonpsych')
print(n_inGnomAD_nonpsych)

print('have_gene')
print(n_have_gene)

print('consequence_category')
print(n_consequence_category)

print('most_severe_consequence')
print(n_most_severe_consequence)

# Note that this is not restricted to the prefiltered variants.
# Let's check the counts filtered to the target intervals:
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(locus=ht_initial_variants.locus, alleles=ht_initial_variants.alleles)
ht = ht.filter(hl.is_defined(ht_initial_variants[ht.key]))

n = ht.count()
n_inGnomAD_nonpsych = ht.aggregate(hl.agg.counter(ht.inGnomAD_nonpsych))
n_have_gene = ht.aggregate(hl.agg.counter(hl.is_defined(ht.vep.worst_csq_for_variant_canonical.gene_symbol)))
n_consequence_category = ht.aggregate(hl.agg.counter(ht.consequence_category))
n_most_severe_consequence = ht.aggregate(hl.agg.counter(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence))

print('n variants:')
print(n)

print('gnomAD_nonpsych')
print(n_inGnomAD_nonpsych)

print('have_gene')
print(n_have_gene)

print('consequence_category')
print(n_consequence_category)

print('most_severe_consequence')
print(n_most_severe_consequence)
