import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

BAM_METRICS = 'gs://dalio_bipolar_w1_w2/data/samples/bam_metrics.tsv'
PHENOFILE = 'gs://dalio_bipolar_w1_w2/data/samples/Dalio_phenotypes.tsv'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes_old_callset.ht'

ht_bam = hl.import_table(BAM_METRICS, key='Samples',
	types={'pct_chimeras': hl.tfloat64, 'pct_contamination': hl.tfloat64})

ht_pheno = hl.import_table(PHENOFILE, key = 'Samples', quote="\"")
ht_pheno = ht_pheno.annotate(Phenotype = ht_pheno.Primary_Disease_Coarse)
ht_pheno = ht_pheno.select('Project', 'Sex', 'PI', 'Phenotype', 'Primary_Disease', 'Location')
ht_pheno = ht_pheno.join(ht_bam, how='left')

ht_pheno.write(PHENOTYPES_TABLE, overwrite=True)

n = ht_pheno.count()

print('')
print('nSamples: ', '{:,}'.format(n))
pprint(ht_pheno.describe())
