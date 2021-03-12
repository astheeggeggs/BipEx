import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict_updated_phenotypes.mt'
GWAS_EXOME_GWAS_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.ht'
GWAS_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.tsv.bgz'
PHENOTYPE_TABLE_BOOL = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_phenotype_table_bool.ht'
GNOMAD_SITES_NON_PSYCH_38_HT = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht'

# Prepare for liftover
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

# BSC results
# If this has not yet been lifted over, need to do that.
# Ensure that it is present in the bucket
# Local location is 
BSC_COUNTS = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/BSC_MAC5_counts.tsv'

def flip_base(base: str) -> str:
    """
    Returns the complement of a base
    :param str base: Base to be flipped
    :return: Complement of input base
    :rtype: str
    """
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))

ht_bsc = hl.import_table(BSC_COUNTS, impute=True)
ht_bsc = ht_bsc.annotate(locus = hl.locus(contig = ht_bsc.chrom, pos = ht_bsc.pos),
    alleles = [ht_bsc.ref, ht_bsc.alt])
ht_bsc = ht_bsc.key_by(ht_bsc.locus, ht_bsc.alleles)

ht_bsc = ht_bsc.annotate(new_locus=hl.liftover(ht_bsc.locus, 'GRCh38', include_strand=True))
ht_bsc = ht_bsc.filter(hl.is_defined(ht_bsc.new_locus))

ht_bsc = ht_bsc.annotate(
        new_alleles = hl.cond(ht_bsc.new_locus.is_negative_strand,
        	[flip_base(ht_bsc.alleles[0]), flip_base(ht_bsc.alleles[1])], ht_bsc.alleles)
        )

ht_bsc = ht_bsc.key_by(locus=ht_bsc.new_locus.result, alleles=ht_bsc.new_alleles)

# Write the result to file.
ht_bsc.write('gs://raw_data_bipolar_dalio_w1_w2/inputs/BSC_MAC5_counts_GRCh38.ht', overwrite=True)

mt = hl.read_matrix_table(QC_MT)

# Create Boolean annotations for all the regressions.
mt = mt.annotate_cols(
    is_BP = hl.case().when(mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder", True)
                    .when(mt.phenotype.PHENOTYPE_COARSE == "Control", False)
                    .default(hl.null(hl.tbool))
    )

mt_bsc = mt.filter_rows(hl.is_defined(ht_bsc[mt.row_key]))
mt_bsc = mt_bsc.repartition(64)

# Filter to the collection of variants following liftover - are there any that are present in both?
# If so, we need to update their allele count.
mt_bsc = mt_bsc.annotate_rows(
    MAC = hl.agg.sum(mt_bsc.GT.n_alt_alleles() > 0),
    burden = hl.agg.sum(mt_bsc.GT.n_alt_alleles())
    )

mt_bsc_bipex = mt_bsc.rows()
mt_bsc_bipex.flatten().export('gs://raw_data_bipolar_dalio_w1_w2/analysis/03_bsc_bipex_intersection_variants.tsv')

# Then, we also want the answer as to whether the variant is in gnomAD.
gnomad_nonpsych_variants_ht = hl.read_table(GNOMAD_SITES_NON_PSYCH_38_HT)
ht_bsc = ht_bsc.annotate(inGnomAD_nonpsych = hl.is_defined(gnomad_nonpsych_variants_ht[ht_bsc.key]))
ht_bsc.flatten().export('gs://raw_data_bipolar_dalio_w1_w2/analysis/03_BSC_MAC5_GRCh38_and_gnomAD.tsv')


