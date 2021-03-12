import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Now we need --requester-pays-allow-all --vep GRCh37 when starting the cluster
# hailctl dataproc start dp --vep GRCh37 --requester-pays-allow-all --region us-central1

# Ensure that the variant list is moved to the cloud
# gsutil cp ~/Repositories/BipEx/BSC_variant_data/BSC_BipEx_gene_variants.v2.txt gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/
# gsutil cp ~/Repositories/BipEx/BSC_data/variant_data/BSC_BipEx_gene_variants.v3.20201018.txt gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/ <- New file with updated definition based on SCHEMA thresholds, MAC <= 5 not in gnomAD top 10.

# Read in the file, and annotate using VEP...making sure that the reference file is build 37.

ht = hl.import_table(#'gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/BSC_BipEx_gene_variants.v2.txt',
	'gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/BSC_BipEx_gene_variants.v3.20201018.txt', 
	impute=True, no_header=True)
ht = ht.rename({'f0':'gene', 'f2':'position', 'f3': 'ref', 'f4': 'alt'})
ht = ht.annotate(chr = ht.f1.replace('chr', ''))
ht = ht.annotate(locus=hl.locus(ht.chr, ht.position, reference_genome = 'GRCh37'),
	alleles = [ht.ref, ht.alt])
ht = ht.key_by(ht.locus, ht.alleles)

# Now create locus and alleles
# Create a compound row key to allow us to annotate.

# Note that the location of the loftee file changed
# (check https://hail.is/docs/0.2/methods/genetics.html?highlight=vep#hail.methods.vep for most up to date location).
ht_vep = hl.vep(ht, "gs://hail-us-vep/vep85-loftee-gcloud.json")
ht_vep.write("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/bsc_variants_vep_annotate.ht", overwrite=True)

