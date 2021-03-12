import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Now we need --requester-pays-allow-all --vep GRCh37 when starting the cluster
# hailctl dataproc start dp --vep GRCh38 --requester-pays-allow-all --region us-central1-b

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'

ht = hl.read_matrix_table(MT).rows()

# vep GRCh38 has moved to gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json
# This is the old version:
# ht_vep = hl.vep(ht, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")
ht_vep = hl.vep(ht, "gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json")
ht_vep.write("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/vep_annotate.ht", overwrite=True)
