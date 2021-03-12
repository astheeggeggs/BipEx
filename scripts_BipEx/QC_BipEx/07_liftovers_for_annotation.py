import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Prepare for liftover
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

# MPC score.
# If this has not yet been lifted over, need to do that.
MPC_SCORE = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2.txt.bgz'

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

mpc_ht = hl.import_table(MPC_SCORE, impute=True)
mpc_ht = mpc_ht.annotate(locus = hl.locus(contig = mpc_ht.chrom, pos = mpc_ht.pos),
    alleles = [mpc_ht.ref, mpc_ht.alt])
mpc_ht = mpc_ht.key_by(mpc_ht.locus, mpc_ht.alleles).select('MPC')

mpc_ht = mpc_ht.annotate(new_locus=hl.liftover(mpc_ht.locus, 'GRCh38', include_strand=True))
mpc_ht = mpc_ht.filter(hl.is_defined(mpc_ht.new_locus))

mpc_ht = mpc_ht.annotate(
        new_alleles = hl.cond(mpc_ht.new_locus.is_negative_strand,
        	[flip_base(mpc_ht.alleles[0]), flip_base(mpc_ht.alleles[1])], mpc_ht.alleles)
        )

mpc_ht = mpc_ht.key_by(locus=mpc_ht.new_locus.result, alleles=mpc_ht.new_alleles)

# Write the result to file.
mpc_ht.write('gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2_GRCh38.ht', overwrite=True)

# Next, need to liftover the entirety of the gnomad data
GNOMAD_SITES_HT = 'gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht'

gnomad_ht = hl.read_table(GNOMAD_SITES_HT)
non_neuro_location = gnomad_ht.freq_index_dict.collect()[0]['non_neuro']
print(non_neuro_location)
gnomad_nonpsych_variants_ht = gnomad_ht.filter(gnomad_ht.freq[non_neuro_location].AF > 0)

# This needs to be lifted over the 38 as well.

gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.annotate(new_locus=hl.liftover(gnomad_nonpsych_variants_ht.locus, 'GRCh38', include_strand=True))
gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.filter(hl.is_defined(gnomad_nonpsych_variants_ht.new_locus))

gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.annotate(
        new_alleles = hl.cond(gnomad_nonpsych_variants_ht.new_locus.is_negative_strand,
        	[flip_base(gnomad_nonpsych_variants_ht.alleles[0]), flip_base(gnomad_nonpsych_variants_ht.alleles[1])], gnomad_nonpsych_variants_ht.alleles)
        )

gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.key_by(locus=gnomad_nonpsych_variants_ht.new_locus.result, alleles=gnomad_nonpsych_variants_ht.new_alleles)
gnomad_nonpsych_variants_ht.write('gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht', overwrite=True)
