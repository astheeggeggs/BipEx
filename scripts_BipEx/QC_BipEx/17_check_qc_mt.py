import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_HARDCALLS_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls.mt'

mt = hl.read_matrix_table(QC_HARDCALLS_MT)

n = mt.count()

print('')
print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

n = mt.filter_rows( (mt.qc.AC[0] == 1) | (mt.qc.AC[1] == 1) ).count()

print('')
print('n singletons:')
print(n[0])

n = mt.filter_rows( (mt.qc.AC[0] == 2) | (mt.qc.AC[1] == 2) ).count()

print('')
print('n doubletons:')
print(n[0])

# This isn't working anymore...weird.
mt = mt.annotate_rows(type = hl.case().when((mt.alleles[0].length() == 1) | (mt.alleles[1].length() == 1), "SNP")
	.when(mt.alleles[0].length() < mt.alleles[1].length(), "Insertion")
	.when(mt.alleles[0].length() > mt.alleles[1].length(), "Deletion")
	.default("No type"))


n = mt.filter_rows(mt.type == "SNP").count()

print('')
print('n SNPS:')
print(n[0])

n = mt.filter_rows(mt.type == "Insertion").count()

print('')
print('n insertions:')
print(n[0])


n = mt.filter_rows(mt.type == "Deletion").count()

print('')
print('n deletions:')
print(n[0])

n = mt.filter_cols(mt.imputesex.impute_sex.is_female).count()

print('')
print('n females:')
print(n[1])

n = mt.filter_cols(~mt.imputesex.impute_sex.is_female).count()

print('')
print('n males:')
print(n[1])

n = mt.filter_cols(mt.phenotype.PHENOTYPE_COARSE == "Control").count()

print('')
print('n controls:')
print(n[1])

n = mt.filter_cols(mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder").count()

print('')
print('n cases:')
print(n[1])

pprint(mt.describe())
