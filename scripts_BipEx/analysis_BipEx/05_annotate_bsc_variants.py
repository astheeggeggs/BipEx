import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# When initialising the cluster, ensure that you use this flag:
# --init gs://gnomad-public/tools/inits/master-init.sh <- this is not used anymore.

# There's a change, now need to include the flag --packages gnomad on cluster start up.
# Also need to add requester pays and update the gnomad public bucket location!
# hailctl dataproc start dp --requester-pays-allow-all --region us-central1 --packages gnomad

from gnomad.utils.vep import process_consequences

# GNOMAD_SITES_HT = 'gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht'
GNOMAD_SITES_HT = 'gs://gnomad-public-requester-pays/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht'
MPC_SCORE = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2.txt.bgz'
GNOMAD_SITES_NON_PSYCH_37_HT = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh37.ht'
# ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/bsc_gene.ht'
ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/bsc_gene_v2.ht'

# Want the non-psych variants list.
gnomad_ht = hl.read_table(GNOMAD_SITES_HT)
non_neuro_location = gnomad_ht.freq_index_dict.collect()[0]['non_neuro']
gnomad_nonpsych_variants_ht = gnomad_ht.filter(gnomad_ht.freq[non_neuro_location].AF > 0)
gnomad_nonpsych_variants_ht.write(GNOMAD_SITES_NON_PSYCH_37_HT, overwrite=True)

# MPC score.
mpc_ht = hl.import_table(MPC_SCORE, impute=True)
mpc_ht = mpc_ht.annotate(locus = hl.locus(contig = mpc_ht.chrom, pos = mpc_ht.pos),
    alleles = [mpc_ht.ref, mpc_ht.alt])
mpc_ht = mpc_ht.key_by(mpc_ht.locus, mpc_ht.alleles).select('MPC')
mpc_ht.write('gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2_GRCh37.ht', overwrite=True)
mpc_ht = hl.read_table('gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2_GRCh37.ht')
mpc_ht = mpc_ht.select(mpc_ht.MPC)

# CADD
# cadd_ht = hl.read_table("gs://hail-datasets-hail-data/CADD.1.4.GRCh37.ht")
cadd_ht = hl.experimental.load_dataset(name='CADD', version='1.4', reference_genome='GRCh37', region='us', cloud='gcp')

# Annotate with the vep information
# ht = hl.read_table("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/bsc_variants_vep_annotate.ht")
ht = hl.read_table("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/bsc_variants_vep_annotate_v2.ht")
ht = ht.annotate(mpc = mpc_ht[ht.key])
ht = ht.annotate(cadd = cadd_ht[ht.key])

ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])

missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])

synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])

non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])

# This time, we want to always use the canonical transcript.

# Here, use Konrad's function from the collection of gnomad functions to annotate get a 
# one-one correspondance between variant and gene/consequence.
ht = process_consequences(ht)

# Find out what the annotation is that I require.
# I will use worst_csq_for_variant_canonical.

ht = ht.annotate(consequence_category = 
    hl.case().when(ptv.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence), "ptv")
             .when(missense.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (~hl.is_defined(ht.vep.worst_csq_for_variant_canonical.polyphen_prediction) | 
                    ~hl.is_defined(ht.vep.worst_csq_for_variant_canonical.sift_prediction) ), "other_missense")
             .when(missense.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (ht.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") & 
                   (ht.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
             .when(missense.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense")
             .when(synonymous.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence), "synonymous")
             .when(non_coding.contains(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence), "non_coding")
             .default("NA")
    )


# Now that we have the gene information for each variant, we can evaluate the constraint.

gnomad_nonpsych_variants_ht = hl.read_table(GNOMAD_SITES_NON_PSYCH_37_HT)
ht = ht.annotate(inGnomAD_nonpsych = hl.is_defined(gnomad_nonpsych_variants_ht[ht.key]))

ht.write(ANNOTATION_TABLE, overwrite=True)

# Now, need to create the annotations required by the collaborators.
# We want 'ptv', 'other_missense', etc 
# inGnomAD_nonpsych
# MPC - throw in for good measure.
ht = hl.read_table(ANNOTATION_TABLE)
ht = ht.key_by().select('gene', 'f1', 'position', 'ref', 'alt', 'consequence_category', 'inGnomAD_nonpsych').rename({'f1': 'chromosome'})
# ht.export('gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/BSC_BipEx_gene_variants.v2.annotated.tsv.bgz')
ht.export('gs://dalio_bipolar_w1_w2_hail_02/data/BSC_variants/BSC_BipEx_gene_variants.v3.20201018.annotated.tsv.bgz')

# Write to disk in exactly the same format that they sent to me (just with a few extra columns).
