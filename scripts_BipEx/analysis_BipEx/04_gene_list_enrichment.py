import hail as hl

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes_rekeyed.mt'

# The following are to generate .tsv file that are gene x sample. These are then read by R and spliced and diced to 
# look for gene set enrichment.

# Using Gene symbol
## TSV files
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC5_gene_counts_per_sample.tsv.bgz'
GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC5_gene_counts_per_sample.tsv.bgz'
# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz'
GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC_gnom_non_psych_gene_counts_per_sample.tsv.bgz'

# Using Gene ID.
## TSV files
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC5_gene_counts_per_sample_ENSG.tsv.bgz'
GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC5_gene_counts_per_sample_ENSG.tsv.bgz'
# Not in GnomAD non psych
GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_BP_including_BPSCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz'
GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/04_SCZ_MAC_gnom_non_psych_gene_counts_per_sample_ENSG.tsv.bgz'

def gene_burden_annotations_per_sample(mt):

    mt = mt.group_rows_by(
        gene_symbol = mt.annotation.vep.worst_csq_for_variant_canonical.gene_symbol,
        consequence_category = mt.annotation.consequence_category
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))

    return mt

def gene_burden_ENSG_annotations_per_sample(mt):

    mt = mt.group_rows_by(
        gene_id = mt.annotation.vep.worst_csq_for_variant_canonical.gene_id,
        consequence_category = mt.annotation.consequence_category
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))

    return mt

mt = hl.read_matrix_table(QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT)
mt = mt.annotate_cols(forest_location = mt['phenotype']['LOCATION'].replace(".*, ", ""))
mt = mt.annotate_cols(forest_location= hl.case()
    .when(((mt['phenotype']['LOCATION'].contains('UK')) | (mt['phenotype']['LOCATION'].contains('IRE'))), "UK/Ireland")
    .when(mt['phenotype']['LOCATION'].contains('Umea, SWE'), "SWE, Umea")
    .when(mt['phenotype']['LOCATION'].contains('Stockholm, SWE'), "SWE, Stockholm")
    .default(mt['forest_location']))

mt_MAC5 = mt.filter_rows(mt.is_MAC5)

# MAC >= 5. Always define MAC based on the entirety of the dataset.

## Bipolar Disorder
bool_phenotypes = ['is_BP1', 'is_BP2', 'is_BP', 'is_BP_including_BPSCZ', 'is_BPNOS', 'is_BPSCZ', 'is_BPPSY', 'is_BP_no_PSY']
mt_BP = mt_MAC5.filter_cols(hl.is_defined(mt_MAC5['phenotype_boolean']['is_BP_including_BPSCZ']))
# mt_gene_BP = gene_burden_annotations_per_sample(mt_BP)
# mt_gene_BP.n.export(GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_TSV)
mt_gene_BP = gene_burden_ENSG_annotations_per_sample(mt_BP)
mt_gene_BP.n.export(GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV)

# Further filter the variants to those not in GnomAD nonpsych.
mt_BP = mt_BP.filter_rows(~mt_BP.annotation.inGnomAD_nonpsych)
# mt_gene_BP = gene_burden_annotations_per_sample(mt_BP)
# mt_gene_BP.n.export(GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV)
mt_gene_BP = gene_burden_ENSG_annotations_per_sample(mt_BP)
mt_gene_BP.n.export(GENE_OUT_BP_INCLUDING_BPSCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV)

## Schizophrenia
bool_phenotypes = ['is_SCZ']
mt_SCZ = mt_MAC5.filter_cols(hl.is_defined(mt_MAC5['phenotype_boolean']['is_SCZ']) &
    ((mt_MAC5.phenotype.LOCATION.contains('UK')) | (mt_MAC5.phenotype.LOCATION.contains('IRE')))
)
# mt_gene_SCZ = gene_burden_annotations_per_sample(mt_SCZ)
# mt_gene_SCZ.n.export(GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_TSV)
mt_gene_SCZ = gene_burden_ENSG_annotations_per_sample(mt_SCZ)
mt_gene_SCZ.n.export(GENE_OUT_SCZ_SAMPLE_MAC5_COUNTS_ENSG_TSV)

# Further filter the variants to those not in GnomAD nonpsych.
mt_SCZ = mt_SCZ.filter_rows(~mt_SCZ.annotation.inGnomAD_nonpsych)
# mt_gene_SCZ = gene_burden_annotations_per_sample(mt_SCZ)
# mt_gene_SCZ.n.export(GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_TSV)
mt_gene_SCZ = gene_burden_ENSG_annotations_per_sample(mt_SCZ)
mt_gene_SCZ.n.export(GENE_OUT_SCZ_SAMPLE_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV)
