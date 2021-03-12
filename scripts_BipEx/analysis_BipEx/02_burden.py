import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_HARDCALLS_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes.mt'
PHENOTYPE_TABLE_BOOL = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_phenotype_table_bool.ht'

# Output (MAC 5).
SAMPLE_MAC5_BURDEN_FILE_BP_including_BPSCZ = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_MAC5_burden_BP_including_BPSCZ.tsv.bgz'
SAMPLE_MAC5_BURDEN_FILE_SCZ = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_MAC5_burden_SCZ.tsv.bgz'

SAMPLE_MAC5_BURDEN_FILE_BP_including_BPSCZ_HT = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_MAC5_burden_BP_including_BPSCZ.ht'
SAMPLE_MAC5_BURDEN_FILE_SCZ_HT = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/02_sample_MAC5_burden_SCZ.ht'

# Fist annotate with the gene constraint information from GnomAD.
mt = hl.read_matrix_table(QC_HARDCALLS_MT)

GENE_CONSTRAINT = "gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
gnom = hl.import_table(GENE_CONSTRAINT, impute=True)
gnom = gnom.key_by(gnom.gene_id)
mt = mt.annotate_rows(gene = gnom[mt.annotation.vep.worst_csq_for_variant_canonical.gene_id])
# Also include the ensembl ID for downstream matching of start and end positions of the gene for manhattan plotting.

# Annotate with the boolean phenotype information.
ht = hl.read_table(PHENOTYPE_TABLE_BOOL)
mt = mt.annotate_cols(phenotype_boolean = ht[mt.col_key])

# Annotate variants with counts from non-psychiatric version of ExAC v0.3.
# Fill in variants not in ExAC variant list with nExAC = 0.
# Annotate variants with LoF/damaging missense annotation.

def burden_annotations(mt, root_ann='burden', annotate=True):
    
    mt = mt.annotate_cols(**{root_ann:hl.struct()})
    ann_data = mt[root_ann].annotate(
        # URVs (Singletons)
        n_URV = hl.agg.count_where(mt.GT.is_non_ref()),
        n_URV_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),
        # Coding Singletons
        n_coding_URV = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category != "non_coding")),
        n_coding_URV_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.annotation.consequence_category != "non_coding")),
        n_coding_URV_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.annotation.consequence_category != "non_coding")),
        # PTVs
        n_URV_PTV = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "ptv")),
        n_URV_PTV_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "ptv") & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_PTV_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "ptv") & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),
        # Damaging missense
        n_URV_damaging_missense = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "damaging_missense")),
        n_URV_damaging_missense_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "damaging_missense") & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_damaging_missense_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "damaging_missense") & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),
        # Other missense
        n_URV_other_missense = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "other_missense")),
        n_URV_other_missense_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "other_missense") & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_other_missense_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "other_missense") & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),
        # Synonymous
        n_URV_synonymous = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "synonymous")),
        n_URV_synonymous_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "synonymous") & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_synonymous_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "synonymous") & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),
        # Non-coding
        n_URV_non_coding = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "non_coding")),
        n_URV_non_coding_SNP = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "non_coding") & (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_non_coding_indel = hl.agg.count_where((mt.GT.is_non_ref()) & (mt.annotation.consequence_category == "non_coding") & (hl.is_indel(mt.alleles[0], mt.alleles[1]))),

        # Next, determine counts with MPC >= 2
        n_URV_MPC_2_damaging_missense = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            (mt.annotation.consequence_category == "damaging_missense") &
            (mt.annotation.mpc.MPC >= 2)),
        n_URV_MPC_2_damaging_missense_SNP = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            (mt.annotation.consequence_category == "damaging_missense") &
            (mt.annotation.mpc.MPC >= 2) &
            (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_MPC_2_damaging_missense_indel = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            (mt.annotation.consequence_category == "damaging_missense") &
            (mt.annotation.mpc.MPC >= 2) &
            (hl.is_indel(mt.alleles[0], mt.alleles[1]))),

        n_URV_MPC_2_missense = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            ((mt.annotation.consequence_category == "damaging_missense") | (mt.annotation.consequence_category == "other_missense"))
            & (mt.annotation.mpc.MPC >= 2)),
        n_URV_MPC_2_missense_SNP = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            ((mt.annotation.consequence_category == "damaging_missense") | (mt.annotation.consequence_category == "other_missense")) &
            (mt.annotation.mpc.MPC >= 2) &
            (hl.is_snp(mt.alleles[0], mt.alleles[1]))),
        n_URV_MPC_2_missense_indel = hl.agg.count_where(
            (mt.GT.is_non_ref()) &
            ((mt.annotation.consequence_category == "damaging_missense") | (mt.annotation.consequence_category == "other_missense")) &
            (mt.annotation.mpc.MPC >= 2) &
            (hl.is_indel(mt.alleles[0], mt.alleles[1])))
    )
    
    if annotate:
        return mt.annotate_cols(**{root_ann: ann_data})
    else:
        return ann_data

# Function to perform the same annotations across the different pLI thresholds.
def burden_across_thresholds(mt, bool_phenotype):
    
    mt = mt.filter_cols(hl.is_defined(mt['phenotype_boolean'][bool_phenotype]))

    mt = burden_annotations(mt)

    mt = mt.filter_rows(~mt.annotation.inGnomAD_nonpsych)
    mt = burden_annotations(mt, root_ann='burden_gnom_non_psych')

    mt_01 = mt.filter_rows(mt.gene.pLI < 0.1)
    mt_01 = burden_annotations(mt_01, 'burden_gnom_non_psych_pli_01')
    mt = mt.annotate_cols(**{'burden_gnom_non_psych_pli_01':hl.struct()})
    mt.annotate_cols(**{'burden_gnom_non_psych_pli_01':mt_01.index_cols(mt.s).burden_gnom_non_psych_pli_01})

    mt_01_09 = mt.filter_rows((mt.gene.pLI > 0.1) & (mt.gene.pLI < 0.9))
    mt_01_09 = burden_annotations(mt_01_09, 'burden_gnom_non_psych_pli_01_09')
    mt = mt.annotate_cols(**{'burden_gnom_non_psych_pli_01_09':hl.struct()})
    mt.annotate_cols(**{'burden_gnom_non_psych_pli_01_09':mt_01_09.index_cols(mt.s).burden_gnom_non_psych_pli_01_09})

    mt = mt.filter_rows(mt.gene.pLI > 0.9)
    mt = burden_annotations(mt, root_ann='burden_gnom_non_psych_pli_09')

    mt = mt.filter_rows(mt.gene.pLI > 0.995)
    mt = burden_annotations(mt, root_ann='burden_gnom_non_psych_pli_0995')

    ht = mt.cols().select(
        'phenotype', 'phenotype_boolean', 'imputesex', 'sample_qc', 'pca',
        'burden', 'burden_gnom_non_psych', 'burden_gnom_non_psych_pli_01',
        'burden_gnom_non_psych_pli_01_09', 'burden_gnom_non_psych_pli_09',
        'burden_gnom_non_psych_pli_0995'
        )
    
    ht = ht.annotate(pca=ht.pca.drop('scores'))

    return(ht)

mt = mt.annotate_rows(
    is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1,
    is_MAC5 = (hl.agg.sum(mt.GT.n_alt_alleles()) <= 5) & (hl.agg.sum(mt.GT.n_alt_alleles()) >= 1)
    )

# MAC5 burden.
mt = mt.filter_rows(mt.is_MAC5)

ht_BP_including_BPSCZ = burden_across_thresholds(mt, 'is_BP_including_BPSCZ')
ht_SCZ = burden_across_thresholds(mt, 'is_SCZ')

ht_BP_including_BPSCZ.write(SAMPLE_MAC5_BURDEN_FILE_BP_including_BPSCZ_HT, overwrite=True)
ht_SCZ.write(SAMPLE_MAC5_BURDEN_FILE_SCZ_HT, overwrite=True)

ht_BP_including_BPSCZ = hl.read_table(SAMPLE_MAC5_BURDEN_FILE_BP_including_BPSCZ_HT)
ht_SCZ = hl.read_table(SAMPLE_MAC5_BURDEN_FILE_SCZ_HT)

ht_BP_including_BPSCZ.flatten().export(SAMPLE_MAC5_BURDEN_FILE_BP_including_BPSCZ)
ht_SCZ.flatten().export(SAMPLE_MAC5_BURDEN_FILE_SCZ)
