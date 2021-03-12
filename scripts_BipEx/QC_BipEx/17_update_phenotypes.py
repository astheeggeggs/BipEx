import hail as hl

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.mt'
QC_HARDCALLS_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls.mt'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

QC_MT_UPDATED = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict_updated_phenotypes.mt'
QC_HARDCALLS_MT_UPDATED = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes.mt'

FINAL_VARIANT_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/17_final_qc.variants.tsv.bgz'
FINAL_SAMPLE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/17_final_qc.samples.tsv.bgz'

ht = hl.read_table(PHENOTYPES_TABLE)

mt = hl.read_matrix_table(QC_HARDCALLS_MT)
mt = mt.drop(mt.phenotype)
mt = mt.annotate_cols(phenotype = ht[mt.s])
mt.write(QC_HARDCALLS_MT_UPDATED, overwrite=True)

mt = hl.read_matrix_table(QC_MT)
mt = mt.drop(mt.phenotype)
mt = mt.annotate_cols(phenotype = ht[mt.s])
mt.write(QC_MT_UPDATED, overwrite=True)

mt = hl.read_matrix_table(QC_HARDCALLS_MT_UPDATED)
ht = mt.cols()

ht.select(
	PROJECT_OR_COHORT = ht.phenotype.PROJECT_OR_COHORT,
	LOCATION = ht.phenotype.LOCATION,
	INSTITUTION = ht.phenotype.INSTITUTION,
	PI = ht.phenotype.PI,
	PHENOTYPE_COARSE = ht.phenotype.PHENOTYPE_COARSE,
	PHENOTYPE_FINE = ht.phenotype.PHENOTYPE_FINE,
	PSYCHOSIS = ht.phenotype.PSYCHOSIS,
	AGE_FI_24 = ht.phenotype.AGE_FI_24,
	AGE_FI_40 = ht.phenotype.AGE_FI_40,
	AGE_FS_24 = ht.phenotype.AGE_FS_24,
	AGE_FS_40 = ht.phenotype.AGE_FS_40,
	AGE_D_24 = ht.phenotype.AGE_D_24,
	AGE_D_40 = ht.phenotype.AGE_D_40,
	PCT_CHIMERAS = ht.phenotype.PCT_CHIMERAS,
	PCT_CONTAMINATION = ht.phenotype.PCT_CONTAMINATION,
	IS_FEMALE = ht.imputesex.impute_sex.is_female,
	GENDER = ht.phenotype.GENDER,
	QC = ht.sample_qc,
	PCA = ht.pca).flatten().export(FINAL_SAMPLE_QC_FILE)

# Want to write out the alternative allele count!
ht = mt.rows()
ht.select(
	rsid = ht.annotation.rsid,
	was_split = ht.annotation.was_split,
	gene_symbol = ht.annotation.vep.worst_csq_for_variant_canonical.gene_symbol,
	gene_id = ht.annotation.vep.worst_csq_for_variant_canonical.gene_id,

	# Still need to add these annotations...latest hg38 for GnomAD is available, can use that I think. Yep.
	# I currently just annotate this after.
	# pli_nopsych = ht.annotation.gene_constraint.pli_nopsych,
	# lof_z_nopsych = ht.annotation.gene_constraint.lof_z_nopsych,
	# mis_z_nopsych = ht.annotation.gene_constraint.mis_z_nopsych,
	# syn_z_nopsych = ht.annotation.gene_constraint.syn_z_nopsych,

	most_severe_consequence = ht.annotation.vep.worst_csq_for_variant_canonical.most_severe_consequence,
	consequence_category = ht.annotation.consequence_category,

	polyphen_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.polyphen_prediction,
	sift_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.sift_prediction,

	inGnomAD_nonpsych = ht.annotation.inGnomAD_nonpsych,
	cadd_phred_score = ht.annotation.cadd.PHRED_score,
	mpc_score = ht.annotation.mpc.MPC,

	qc = ht.qc).flatten().export(FINAL_VARIANT_QC_FILE)
 