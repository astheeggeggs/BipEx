import hail as hl

from hail.plot import show
from pprint import pprint
import numpy as np
import pandas as pd

hl.plot.output_notebook()

QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes_rekeyed.mt'
GENE_SETS_BP_INCLUDING_BPSCZ_MAC5_GNOM_NON_PSYCH_COUNTS_TSV = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.tsv"

def run_regressions_obs(mt, phenotypes):

	init = True
	for phenotype in phenotypes:

		logreg_ht = hl.logistic_regression_rows(test = 'lrt',
			y = mt.phenotype_boolean[phenotype],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10],
			pass_through = ['gene_set', 'consequence_category']
		).annotate(phenotype=phenotype)

		firth_ht = hl.logistic_regression_rows(test = 'firth',
			y = mt.phenotype_boolean[phenotype],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10],
			pass_through = ['gene_set', 'consequence_category']
		).annotate(phenotype=phenotype)
		
		mt = mt.annotate_rows(linreg = hl.agg.linreg(y=mt.x,
			x=[1, mt.phenotype_boolean[phenotype], mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10]))
		linreg_ht = mt.rows().annotate(phenotype=phenotype)
		linreg_ht = linreg_ht.select(
			beta = linreg_ht.linreg.beta[1], 
			standard_error = linreg_ht.linreg.standard_error[1],  
	        t_stat = linreg_ht.linreg.t_stat[1], 
	        p_value = linreg_ht.linreg.p_value[1],
	        phenotype = phenotype)

		if init is True:
			results_logreg_ht = logreg_ht
			results_firth_ht = firth_ht
			results_linreg_ht = linreg_ht
			init = False
		else:
			results_logreg_ht = results_logreg_ht.union(logreg_ht)
			results_firth_ht = results_firth_ht.union(firth_ht)
			results_linreg_ht = results_linreg_ht.union(linreg_ht)

	return results_linreg_ht, results_logreg_ht, results_firth_ht

# For this, need to ensure that coding burden is added as a column (sample) annotation.
def run_regressions_obs_coding(mt, phenotypes):

	init = True
	for phenotype in phenotypes:

		logreg_ht = hl.logistic_regression_rows(test = 'lrt',
			y = mt.phenotype_boolean[phenotype],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.coding_burden],
			pass_through = ['gene_set', 'consequence_category']
		).annotate(phenotype=phenotype)

		firth_ht = hl.logistic_regression_rows(test = 'firth',
			y = mt.phenotype_boolean[phenotype],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.coding_burden],
			pass_through = ['gene_set', 'consequence_category']
		).annotate(phenotype=phenotype)
		
		mt = mt.annotate_rows(linreg = hl.agg.linreg(y=mt.x,
			x=[1, mt.phenotype_boolean[phenotype], mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.coding_burden]))
		linreg_ht = mt.rows().annotate(phenotype=phenotype)
		linreg_ht = linreg_ht.select(
			beta = linreg_ht.linreg.beta[1], 
			standard_error = linreg_ht.linreg.standard_error[1],  
	        t_stat = linreg_ht.linreg.t_stat[1], 
	        p_value = linreg_ht.linreg.p_value[1],
	        phenotype = phenotype)

		if init is True:
			results_logreg_ht = logreg_ht
			results_firth_ht = firth_ht
			results_linreg_ht = linreg_ht
			init = False
		else:
			results_logreg_ht = results_logreg_ht.union(logreg_ht)
			results_firth_ht = results_firth_ht.union(firth_ht)
			results_linreg_ht = results_linreg_ht.union(linreg_ht)

	return results_linreg_ht, results_logreg_ht, results_firth_ht

def convert_gene_set_tsv_to_mt(gene_set_tsv):
	return hl.import_matrix_table(gene_set_tsv,
		row_fields={
		'gene_set': hl.tstr,
		'consequence_category': hl.tstr
		}, delimiter='\t').key_rows_by('gene_set', 'consequence_category')

def prepare_mt(gene_set_tsv, mt_path, mt_outpath):

	mt = convert_gene_set_tsv_to_mt(gene_set_tsv)
	
	ht = hl.read_matrix_table(mt_path).cols().select(
		'phenotype_boolean', 'imputesex', 'sample_qc', 'pca')
	ht = ht.select(
		phenotype_boolean = ht.phenotype_boolean,
		is_female = ht.imputesex.impute_sex.is_female,
		PC1 = ht.pca.PC1,
		PC2 = ht.pca.PC2,
		PC3 = ht.pca.PC3,
		PC4 = ht.pca.PC4,
		PC5 = ht.pca.PC5,
		PC6 = ht.pca.PC6,
		PC7 = ht.pca.PC7,
		PC8 = ht.pca.PC8,
		PC9 = ht.pca.PC9,
		PC10 = ht.pca.PC10)

	mt = mt.annotate_cols(**ht[mt.col_key])
	mt.write(mt_outpath, overwrite=True)

def run_regressions_obs_control_for_burden(mt_full, phenotypes, coding=False):

	if coding is True:
		mt_full = mt_full.filter_rows(mt_full.consequence_category != 'non_coding')

	gene_sets = mt_full.gene_set.collect()
	gene_sets = np.unique(gene_sets)

	init = True

	for gene_set in gene_sets:

		mt = mt_full.filter_rows(mt_full.gene_set == gene_set)
		mt = mt.annotate_cols(burden = hl.agg.sum(mt.x))
		
		for phenotype in phenotypes:

			logreg_ht = hl.logistic_regression_rows(test = 'lrt',
				y = mt.phenotype_boolean[phenotype],
				x = mt.x,
				covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.burden],
				pass_through = ['gene_set', 'consequence_category']
			).annotate(phenotype=phenotype)

			firth_ht = hl.logistic_regression_rows(test = 'firth',
				y = mt.phenotype_boolean[phenotype],
				x = mt.x,
				covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.burden],
				pass_through = ['gene_set', 'consequence_category']
			).annotate(phenotype=phenotype)
			
			mt = mt.annotate_rows(linreg = hl.agg.linreg(y=mt.x, x=[1, mt.phenotype_boolean[phenotype], mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10, mt.burden]))
			linreg_ht = mt.rows().annotate(phenotype=phenotype)
			linreg_ht = linreg_ht.select(
				beta = linreg_ht.linreg.beta[1], 
				standard_error = linreg_ht.linreg.standard_error[1],  
		        t_stat = linreg_ht.linreg.t_stat[1], 
		        p_value = linreg_ht.linreg.p_value[1],
		        phenotype = phenotype)

			if init is True:
				results_logreg_ht = logreg_ht
				results_firth_ht = firth_ht
				results_linreg_ht = linreg_ht
				init = False
			else:
				results_logreg_ht = results_logreg_ht.union(logreg_ht)
				results_firth_ht = results_firth_ht.union(firth_ht)
				results_linreg_ht = results_linreg_ht.union(linreg_ht)

	return results_linreg_ht, results_logreg_ht, results_firth_ht

def permute_phenotypes(np_pheno, n_perms):
	
	np_pheno[np_pheno == None] = 2
	np_pheno = np_pheno.astype(int)
	np_pheno_mat = np.repeat(np_pheno, n_perms).reshape(np_pheno.size, n_perms)
	for i in range(np_pheno_mat.shape[1]):
	    np.random.shuffle(np_pheno_mat[:,i])

	return(np_pheno_mat)


def run_regressions_perm(mt, phenotypes, n_perms):

	mt = mt.add_col_index()

	init = True
	for phenotype in phenotypes:
		
		mt = mt.annotate_globals(
			pheno_perm = permute_phenotypes(np.array(mt.phenotype_boolean[phenotype].collect()), n_perms))

		mt = mt.annotate_globals(
            pheno_perm = mt.pheno_perm.map(
                lambda x: { hl.array([False, True, hl.null(hl.tbool)])[hl.int(x)] }
            )
		)
		
		for perm in range(n_perms):

			logreg_ht = hl.logistic_regression_rows(test = 'lrt',
			y = mt.pheno_perm[mt.col_idx, perm],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3,
				mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10],
			pass_through = ['gene_set', 'consequence_category']
			).annotate(phenotype=phenotype, permutation=perm+1)

			logreg_ht = logreg_ht.select(
				p_value = logreg_ht.p_value,
				beta = logreg_ht.beta,
				chi_sq_stat = logreg_ht.chi_sq_stat,
				permutation = logreg_ht.permutation
			)

			firth_ht = hl.logistic_regression_rows(test = 'firth',
			y = mt.pheno_perm[mt.col_idx, perm],
			x = mt.x,
			covariates = [1, mt.is_female, mt.PC1, mt.PC2, mt.PC3,
				mt.PC4, mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10],
			pass_through = ['gene_set', 'consequence_category']
			).annotate(phenotype=phenotype, permutation=perm+1)
			
			firth_ht = firth_ht.select(
				p_value = logreg_ht.p_value,
				beta = logreg_ht.beta,
				chi_sq_stat = logreg_ht.chi_sq_stat,
				permutation = logreg_ht.permutation
			)

			mt = mt.annotate_rows(linreg = hl.agg.linreg(y=mt.x,
				x=[1, mt.pheno_perm[mt.col_idx, perm], 
				mt.is_female, mt.PC1, mt.PC2, mt.PC3, mt.PC4,
				mt.PC5, mt.PC6, mt.PC7, mt.PC8, mt.PC9, mt.PC10]))

			linreg_ht = mt.rows().annotate(phenotype=phenotype, permutation=perm+1)
			linreg_ht = linreg_ht.select(
				beta = linreg_ht.linreg.beta[1], 
				standard_error = linreg_ht.linreg.standard_error[1],  
				t_stat = linreg_ht.linreg.t_stat[1], 
				p_value = linreg_ht.linreg.p_value[1]
			)

			if init is True:
				perms_logreg_ht = logreg_ht
				perms_firth_ht = firth_ht
				perms_linreg_ht = linreg_ht
				init_pheno = False
			else:
				perms_logreg_ht = perms_logreg_ht.union(logreg_ht)
				perms_firth_ht = perms_firth_ht.union(firth_ht)
				perms_linreg_ht = perms_linreg_ht.union(linreg_ht)

	return perms_linreg_ht, perms_logreg_ht, perms_firth_ht

gene_set_tsv_BP = GENE_SETS_BP_INCLUDING_BPSCZ_MAC5_GNOM_NON_PSYCH_COUNTS_TSV
mt_outpath_BP = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.mt"
folder_BP = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/"
gene_set_obs_output_BP = "BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.ht"
gene_set_obs_output_BP_tsv = "BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists.tsv.bgz"
mt_path = QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT
burden_tsv = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_burden.tsv"
burden_syn_tsv = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_synonymous_burden.tsv"

prepare_mt(gene_set_tsv_BP, mt_path, mt_outpath_BP)

phenotypes_BP = ['is_BP1', 'is_BP2', 'is_BP', 'is_BPPSY', 'is_BP_no_PSY']

mt = hl.read_matrix_table(mt_outpath_BP)
ht_burden = hl.import_table(burden_tsv, impute=True)
ht_burden_syn = hl.import_table(burden_tsv, impute=True)

ht_burden = ht_burden.key_by("sample")
ht_burden_syn = ht_burden_syn.key_by("sample")
mt = mt.annotate_cols(coding_burden=ht_burden[mt.col_key].coding_burden)

linreg_ht, logreg_ht, firth_ht = run_regressions_obs_coding(mt, phenotypes_BP)
linreg_ht.write(folder_BP + 'linreg_coding_burden_' + gene_set_obs_output_BP, overwrite=True)
logreg_ht.write(folder_BP + 'logreg_coding_burden_' + gene_set_obs_output_BP, overwrite=True)
firth_ht.write(folder_BP + 'firth_coding_burden_' + gene_set_obs_output_BP, overwrite=True)

# Now, resave these in tsv files for plotting.
linreg_ht.export(folder_BP + 'linreg_coding_burden_' + gene_set_obs_output_BP_tsv)
logreg_ht.export(folder_BP + 'logreg_coding_burden_' + gene_set_obs_output_BP_tsv)
firth_ht.export(folder_BP + 'firth_coding_burden_' + gene_set_obs_output_BP_tsv)

# And now run the same regressions for the final gene-sets - probing the overlap with SCZ further.

GENE_SETS_BP_INCLUDING_BPSCZ_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_SCZ_overlap_gene_lists.tsv"
gene_set_tsv_BP = GENE_SETS_BP_INCLUDING_BPSCZ_MAC5_GNOM_NON_PSYCH_COUNTS_ENSG_TSV

mt_outpath_BP = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists_SCZ_overlap.mt"
folder_BP = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/output/obs/BP_including_BPSCZ/"
gene_set_obs_output_BP = "BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists_SCZ_overlap.ht"
gene_set_obs_output_BP_tsv = "BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_gene_lists_SCZ_overlap.tsv.bgz"
mt_path = QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT
burden_tsv = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_SCZ_overlap_burden.tsv"
burden_syn_tsv = "gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_sets/BP_including_BPSCZ_MAC5_gnom_non_psych_gene_set_counts_per_sample_SCZ_overlap_synonymous_burden.tsv"

prepare_mt(gene_set_tsv_BP, mt_path, mt_outpath_BP)

phenotypes_BP = ['is_BP1', 'is_BP2', 'is_BP', 'is_BP_including_BPSCZ', 'is_BPNOS', 'is_BPSCZ', 'is_BPPSY', 'is_BP_no_PSY']

mt = hl.read_matrix_table(mt_outpath_BP)
ht_burden = hl.import_table(burden_tsv, impute=True)
ht_burden_syn = hl.import_table(burden_syn_tsv, impute=True)

ht_burden = ht_burden.key_by("sample")
ht_burden_syn = ht_burden_syn.key_by("sample")
mt = mt.annotate_cols(coding_burden=ht_burden[mt.col_key].coding_burden)

linreg_ht, logreg_ht, firth_ht = run_regressions_obs_coding(mt, phenotypes_BP)
linreg_ht.write(folder_BP + 'linreg_coding_burden_' + gene_set_obs_output_BP, overwrite=True)
logreg_ht.write(folder_BP + 'logreg_coding_burden_' + gene_set_obs_output_BP, overwrite=True)
firth_ht.write(folder_BP + 'firth_coding_burden_' + gene_set_obs_output_BP, overwrite=True)

# Now, resave these in tsv files for plotting.
linreg_ht.export(folder_BP + 'linreg_coding_burden_' + gene_set_obs_output_BP_tsv)
logreg_ht.export(folder_BP + 'logreg_coding_burden_' + gene_set_obs_output_BP_tsv)
firth_ht.export(folder_BP + 'firth_coding_burden_' + gene_set_obs_output_BP_tsv)

mt = mt.annotate_cols(coding_burden=ht_burden_syn[mt.col_key].coding_burden)
