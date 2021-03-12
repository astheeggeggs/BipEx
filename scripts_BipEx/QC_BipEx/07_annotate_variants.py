import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# When initialising the cluster, ensure that you use this flag:
# --init gs://gnomad-public/tools/inits/master-init.sh

# There's a change, now need to include the flag --packages gnomad on cluster start up.

# from gnomad_hail import * <- this changed to gnomad.utils.
from gnomad.utils import *

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'

# Want the non-psych variants list.
GNOMAD_SITES_NON_PSYCH_38_HT = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/gnomad.exomes.r2.1.1.non_psych_sites_GRCh38.ht'
ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

# MPC score.
# If this has not yet been lifted over, need to do that.
MPC_SCORE = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2_GRCh38.ht'
mpc_ht = hl.read_table(MPC_SCORE)
mpc_ht = mpc_ht.select(mpc_ht.MPC)

# CADD
cadd_ht = hl.read_table("gs://hail-datasets-hail-data/CADD.1.4.GRCh38.ht")

mt = hl.read_matrix_table(MT)

# Annotate with the vep information
vep_ht = hl.read_table("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/vep_annotate.ht")
vep_ht = vep_ht.select('vep')
mt = mt.annotate_rows(vep_ann = vep_ht[mt.row_key])
mt = mt.annotate_rows(vep = mt.vep_ann.vep)
mt = mt.drop(mt.vep_ann)

mt = mt.annotate_rows(mpc = mpc_ht[mt.row_key])
mt = mt.annotate_rows(cadd = cadd_ht[mt.row_key])

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
mt = process_consequences(mt)

# Find out what the annotation is that I require.
# I will use worst_csq_for_variant_canonical.

mt = mt.annotate_rows(consequence_category = 
    hl.case().when(ptv.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "ptv")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (~hl.is_defined(mt.vep.worst_csq_for_variant_canonical.polyphen_prediction) | 
                    ~hl.is_defined(mt.vep.worst_csq_for_variant_canonical.sift_prediction) ), "other_missense")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (mt.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") & 
                   (mt.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense")
             .when(synonymous.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "synonymous")
             .when(non_coding.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "non_coding")
             .default("NA")
    )


# Now that we have the gene information for each variant, we can evaluate the constraint.

gnomad_nonpsych_variants_ht = hl.read_table(GNOMAD_SITES_NON_PSYCH_38_HT)
mt = mt.annotate_rows(inGnomAD_nonpsych = hl.is_defined(gnomad_nonpsych_variants_ht[mt.row_key]))

mt_rows = mt.rows().repartition(64)
mt_rows.write(ANNOTATION_TABLE, overwrite=True)



