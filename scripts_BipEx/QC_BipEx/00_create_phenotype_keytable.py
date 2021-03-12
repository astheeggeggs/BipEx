import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Older version (which was used to run the QC pipeline).
# PHENOFILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/BIP_phenotype_information_cleaned_new_subtype_information_added.tsv' 
# Note that in the older version I used the 'PHENOTYPE_COARSE_NEW' and 'PHENOTYPE_FINE_NEW' columns. These don't exist in this file (I deleted them to avoid confusion).
# PHENOFILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_final.tsv'
# Updated to correct psychosis and include age of onset information
PHENOFILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_and_psychosis_and_aao_final.tsv'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

pheno_table = hl.import_table(PHENOFILE, quote="\"", key='SAMPLE_ALIAS',
	types={'PCT_CHIMERAS': hl.tfloat64, 'PCT_CONTAMINATION': hl.tfloat64}, impute=True)

pheno_table.select(pheno_table.PROJECT_OR_COHORT, pheno_table.GENDER,
		pheno_table.PI, pheno_table.PHENOTYPE_COARSE, pheno_table.PHENOTYPE_FINE,
		pheno_table.LOCATION, pheno_table.INSTITUTION,
		pheno_table.PCT_CONTAMINATION, pheno_table.PCT_CHIMERAS, pheno_table.PSYCHOSIS,
		pheno_table.AGE_FI_24, pheno_table.AGE_FI_40,
		pheno_table.AGE_FS_24, pheno_table.AGE_FS_40,
		pheno_table.AGE_D_24, pheno_table.AGE_D_40).write(PHENOTYPES_TABLE, overwrite=True)

ht = hl.read_table(PHENOTYPES_TABLE)
n_samples = ht.count()

print('')
print('nSamples: ', '{:,}'.format(n_samples))
pprint(ht.describe())
