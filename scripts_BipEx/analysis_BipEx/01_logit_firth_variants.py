import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict_updated_phenotypes.mt'
GWAS_EXOME_GWAS_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.ht'
GWAS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.tsv.bgz'
PHENOTYPE_TABLE_BOOL = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_phenotype_table_bool.ht'

mt = hl.read_matrix_table(QC_MT)

# Create Boolean annotations for all the regressions.
mt = mt.annotate_cols(
    is_BP1 = hl.case().when(mt.phenotype.PHENOTYPE_FINE == "Bipolar Disorder 1", True)
                    .when(mt.phenotype.PHENOTYPE_FINE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_BP2 = hl.case().when(mt.phenotype.PHENOTYPE_FINE == "Bipolar Disorder 2", True)
                    .when(mt.phenotype.PHENOTYPE_FINE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_BPNOS = hl.case().when(mt.phenotype.PHENOTYPE_FINE == "Bipolar Disorder NOS", True)
                    .when(mt.phenotype.PHENOTYPE_FINE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_BPSCZ = hl.case().when(mt.phenotype.PHENOTYPE_FINE == "Schizoaffective", True)
                    .when(mt.phenotype.PHENOTYPE_FINE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_BP = hl.case().when(mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder", True)
                    .when(mt.phenotype.PHENOTYPE_COARSE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_BP_including_BPSCZ = hl.case().when(((mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder") | (mt.phenotype.PHENOTYPE_COARSE == "Schizoaffective")), True)
                      .when(mt.phenotype.PHENOTYPE_COARSE == "Control", False)
                      .default(hl.null(hl.tbool)),
    is_SCZ = hl.case().when(mt.phenotype.PHENOTYPE_COARSE == "Schizophrenia", True)
                    .when(mt.phenotype.PHENOTYPE_COARSE == "Control", False)
                    .default(hl.null(hl.tbool)),
    is_PSYCHOSIS = mt.phenotype.PSYCHOSIS
    )

mt = mt.annotate_cols(
    is_BPPSY = hl.case().when((mt.is_BP_including_BPSCZ) & (mt.is_PSYCHOSIS), True)
                    .when(~mt.is_BP_including_BPSCZ, False)
                    .default(hl.null(hl.tbool)),
    is_BP_no_PSY = hl.case().when((mt.is_BP_including_BPSCZ) & (~mt.is_PSYCHOSIS), True)
                    .when(~mt.is_BP_including_BPSCZ, False)
                    .default(hl.null(hl.tbool))                
    )

mt.cols().select('is_BP1', 'is_BP2', 'is_BPNOS', 'is_BPSCZ', 'is_BP', 'is_BP_including_BPSCZ',
    'is_SCZ', 'is_BPPSY', 'is_BP_no_PSY', 'is_PSYCHOSIS').write(PHENOTYPE_TABLE_BOOL, overwrite=True)

mt = mt.annotate_rows(MAC = hl.min(hl.agg.sum(mt.GT.n_alt_alleles()), hl.agg.sum(hl.int64(mt.GT.is_het_ref()) + 2 * hl.int64(mt.GT.is_hom_ref()))))
mt_MAC10 = mt.filter_rows(mt.MAC >= 10)

def run_logistic_bool(mt, variable):

    ht = hl.logistic_regression_rows(
        test='firth',
        y = mt[variable], 
        x = mt.GT.n_alt_alleles(),
        covariates=[1, mt.imputesex.impute_sex.is_female,
            mt.pca.PC1, mt.pca.PC2, mt.pca.PC3, mt.pca.PC4, mt.pca.PC5,
            mt.pca.PC6, mt.pca.PC7, mt.pca.PC8, mt.pca.PC9, mt.pca.PC10]
        )

    mt = mt.filter_cols(hl.is_defined(mt[variable]))
    mt = mt.annotate_rows(MAC = hl.min(hl.agg.sum(mt.GT.n_alt_alleles()), hl.agg.sum(hl.int64(mt.GT.is_het_ref()) + 2 * hl.int64(mt.GT.is_hom_ref()))))
    ht = ht.annotate(MAC = mt.rows()[ht.key].MAC)
    return(ht)

ht_BP1 = run_logistic_bool(mt_MAC10, 'is_BP1')
ht_BP2 = run_logistic_bool(mt_MAC10, 'is_BP2')
ht_BPNOS = run_logistic_bool(mt_MAC10, 'is_BPNOS')
ht_BPSCZ = run_logistic_bool(mt_MAC10, 'is_BPSCZ')
ht_BP = run_logistic_bool(mt_MAC10, 'is_BP')
ht_BP_including_BPSCZ = run_logistic_bool(mt_MAC10, 'is_BP_including_BPSCZ')
ht_BPPSY = run_logistic_bool(mt_MAC10, 'is_BPPSY')
ht_BP_no_PSY = run_logistic_bool(mt_MAC10, 'is_BP_no_PSY')

mt_MAC10 = mt_MAC10.annotate_rows(
    logreg = hl.struct(
        BP1 = ht_BP1[mt_MAC10.row_key],
        BP2 = ht_BP2[mt_MAC10.row_key],
        BPNOS = ht_BPNOS[mt_MAC10.row_key],
        BPSCZ = ht_BPSCZ[mt_MAC10.row_key],
        BP = ht_BP[mt_MAC10.row_key],
        BP_including_BPSCZ = ht_BP_including_BPSCZ[mt_MAC10.row_key],
        BPPSY = ht_BPPSY[mt_MAC10.row_key],
        BP_no_PSY = ht_BP_no_PSY[mt_MAC10.row_key]
    ),
    gene_symbol = mt_MAC10.annotation.vep.worst_csq_for_variant_canonical.gene_symbol,
    consequence_category = mt_MAC10.annotation.consequence_category,
    MAF = hl.min(mt_MAC10.qc.AF[0], mt_MAC10.qc.AF[1])
)

mt_MAC10.rows().select(
    'rsid', 'logreg', 'gene_symbol', 'consequence_category', 'MAF', 'MAC'
    ).write(GWAS_EXOME_GWAS_TABLE, overwrite=True)

# Also, write the Bipolar (including Schizoaffective) to two .tsv file for Manhattan plotting.
hl.read_table(GWAS_EXOME_GWAS_TABLE).flatten().export(GWAS_TSV)


