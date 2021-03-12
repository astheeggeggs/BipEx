import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

ds = hl.import_vcf('gs://raw_data_bipolar_dalio_w1_w2_hail_02/dsp_cloud_output/*.hard_filtered_with_genotypes.vcf.gz',
	reference_genome='GRCh38', force_bgz=True, find_replace=('nul', '.'))
ds.write(output='gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes.mt', overwrite=True)
