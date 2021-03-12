import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes_rekeyed.mt'
mt = hl.read_matrix_table(QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT)

# Initially this is just to create the bare minimum of what's required. We will be adding to this as results come in.

# First, select the variant annotations and transform to the required format
BROWSER_VARIANT_ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02_browser/data/ht/browser_variant_annotation_table.ht'
BROWSER_VARIANT_RESULTS_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02_browser/data/ht/browser_variant_results_table.ht'
BROWSER_GENE_RESULTS_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02_browser/data/ht/browser_gene_results_table.ht'

ht = mt.rows().key_by('locus', 'alleles')
# Throw away everything except what we need for variant annotations

ht = ht.select(transcript_csq = ht.annotation.vep.transcript_consequences,
               worst_csq_for_variant_canonical = ht.annotation.vep.worst_csq_for_variant_canonical,
               csq_analysis = ht.annotation.consequence_category,
               csq_worst = ht.annotation.vep.most_severe_consequence,
               cadd = ht.annotation.cadd.PHRED_score,
               mpc = ht.annotation.mpc.MPC)

# Explode on consequence category
ht = ht.explode('transcript_csq', name='transcript_csq')

ht = ht.select(variant_id = hl.variant_str(ht.locus, ht.alleles),
               gene_id = ht.worst_csq_for_variant_canonical.gene_id,
               gene_name = ht.worst_csq_for_variant_canonical.gene_symbol,
               canonical_transcript_id = ht.worst_csq_for_variant_canonical.transcript_id,
               transcript_id = ht.transcript_csq.transcript_id,
               hgvsc_canonical = ht.worst_csq_for_variant_canonical.hgvsc,
               hgvsc = ht.transcript_csq.hgvsc,
               hgvsp_canonical = ht.worst_csq_for_variant_canonical.hgvsp,
               hgvsp = ht.transcript_csq.hgvsp,
               csq_analysis = ht.csq_analysis,
               csq_worst = ht.csq_worst,
               csq_canonical = ht.worst_csq_for_variant_canonical.most_severe_consequence,
               cadd = ht.cadd,
               mpc = ht.mpc,
               polyphen = ht.worst_csq_for_variant_canonical.polyphen_prediction).write(BROWSER_VARIANT_ANNOTATION_TABLE, overwrite=True)

# Next, create the variants results table.
# First, filter down to the set of Bipolar cases and controls (remove the schizophrenia samples).

mt = mt.filter_cols(hl.is_defined(mt.phenotype_boolean.is_BP_including_BPSCZ))
mt = mt.annotate_rows(BP1_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BP1, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BP2_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BP2, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BPNOS_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BPNOS, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BPSCZ_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BPSCZ, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BP_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BP, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BP_including_BPSCZ_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BP_including_BPSCZ, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BPPSY_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BPPSY, hl.agg.call_stats(mt.GT, mt.alleles)),
                      BP_no_PSY_stats = hl.agg.filter(
                        mt.phenotype_boolean.is_BP_no_PSY, hl.agg.call_stats(mt.GT, mt.alleles)),
                      CONTROL_stats = hl.agg.filter(
                        ~mt.phenotype_boolean.is_BP_including_BPSCZ, hl.agg.call_stats(mt.GT, mt.alleles))
                      )

# Read in the variant results and annotate the table with the relevant information
GWAS_EXOME_GWAS_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/BipEx_GWAS.ht'
ht_variant_assoc = hl.read_table(GWAS_EXOME_GWAS_TABLE)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.annotate_rows(**ht_variant_assoc[mt.row_key])

mt = mt.annotate_rows(in_analysis = mt.is_MAC5)
mt = mt.annotate_rows(explode_data = [
        hl.struct(analysis_group='Bipolar Disorder 1',
          case=mt.BP1_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BP1.beta, chi_sq_stat=mt.logreg.BP1.chi_sq_stat,
          p_value=mt.logreg.BP1.p_value),
        hl.struct(analysis_group='Bipolar Disorder 2',
          case=mt.BP2_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BP2.beta, chi_sq_stat=mt.logreg.BP2.chi_sq_stat,
          p_value=mt.logreg.BP2.p_value),
        # hl.struct(analysis_group='Bipolar Disorder NOS',
        #   case=mt.BPNOS_stats, control=mt.CONTROL_stats,
        #   estimate=mt.logreg.BPNOS.beta, chi_sq_stat=mt.logreg.BPNOS.chi_sq_stat,
        #   p_value=mt.logreg.BPNOS.p_value),
        # hl.struct(analysis_group='Schizoaffective',
        #   case=mt.BPSCZ_stats, control=mt.CONTROL_stats,
        #   estimate=mt.logreg.BPSCZ.beta, chi_sq_stat=mt.logreg.BPSCZ.chi_sq_stat,
        #   p_value=mt.logreg.BPSCZ.p_value),
        hl.struct(analysis_group='Bipolar Disorder',
          case=mt.BP_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BP.beta, chi_sq_stat=mt.logreg.BP.chi_sq_stat,
          p_value=mt.logreg.BP.p_value),
        hl.struct(analysis_group='Bipolar Disorder (including Schizoaffective)',
          case=mt.BP_including_BPSCZ_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BP_including_BPSCZ.beta, chi_sq_stat=mt.logreg.BP_including_BPSCZ.chi_sq_stat,
          p_value=mt.logreg.BP_including_BPSCZ.p_value),
        hl.struct(analysis_group='Bipolar Disorder with Psychosis',
          case=mt.BPPSY_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BPPSY.beta, chi_sq_stat=mt.logreg.BPPSY.chi_sq_stat,
          p_value=mt.logreg.BPPSY.p_value),
        hl.struct(analysis_group='Bipolar Disorder without Psychosis',
          case=mt.BP_no_PSY_stats, control=mt.CONTROL_stats,
          estimate=mt.logreg.BP_no_PSY.beta, chi_sq_stat=mt.logreg.BP_no_PSY.chi_sq_stat,
          p_value=mt.logreg.BP_no_PSY.p_value)
        ]
    )

mt = mt.explode_rows('explode_data')
mt = mt.transmute_rows(**mt.explode_data)
ht = mt.rows()
ht = ht.select(
    variant_id = hl.variant_str(ht.locus, ht.alleles),
    analysis_group = ht.analysis_group, 
    ac_case = ht.case.AC,
    an_case = ht.case.AN,
    af_case = ht.case.AF,
    ac_ctrl = ht.control.AC,
    an_ctrl = ht.control.AN,
    af_ctrl = ht.control.AF,
    estimate = ht.estimate,
    chi_sq_stat = ht.chi_sq_stat,
    p_value = ht.p_value,
    in_analysis = ht.in_analysis
)

ht.write(BROWSER_VARIANT_RESULTS_TABLE, overwrite=True)

# GENE_COUNT_BP_including_BPSCZ_MAC5 = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_counts.ht'
GENE_BURDEN_BP_including_BPSCZ_MAC5 = 'gs://dalio_bipolar_w1_w2_hail_02/analysis/03_BP_including_BPSCZ_MAC5_gene_burdens.ht'
# This gives all the counts - need to reannotate this with the p-values from the fisher's exact and CMH tests.

# The last piece is to create the gene results table.
ht = hl.read_table(GENE_BURDEN_BP_including_BPSCZ_MAC5)
ht = ht.filter((ht.consequence_category == "ptv") | (ht.consequence_category == "damaging_missense"))

# Ensure all of the tsv files are moved to the cloud.
# gsutil cp ~/Repositories/BipEx/analysis_plots/gene_counts_qq/plots/BP_*fisher*tsv gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_fisher_tables/
# gsutil cp ~/Repositories/BipEx/analysis_plots/gene_counts_qq/plots/BP_*CMH*tsv gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_CMH_tables/

# First the Fisher's exact for MAC5 (can be in gnomAD)
ht_fisher = hl.import_table("gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_fisher_tables/BP_gene_fisher_MAC5_qq.tsv", impute=True).key_by('gene_symbol', 'consequence_category')
ht_fisher = ht_fisher.select(
  # p-values
  'is_BP1.pval',
  'is_BP2.pval',
  'is_BP.pval',
  'is_BP_including_BPSCZ.pval',
  'is_BPNOS.pval',
  'is_BPSCZ.pval',
  'is_BPPSY.pval',
  'is_BP_no_PSY.pval',
  # Odds ratios
  'is_BP1.OR',
  'is_BP2.OR',
  'is_BP.OR',
  'is_BP_including_BPSCZ.OR',
  'is_BPNOS.OR',
  'is_BPSCZ.OR',
  'is_BPPSY.OR',
  'is_BP_no_PSY.OR'
)

ht_fisher = ht_fisher.rename(
  {
  # p-values
  'is_BP1.pval' : 'is_BP1_pval',
  'is_BP2.pval' : 'is_BP2_pval',
  'is_BP.pval' : 'is_BP_pval',
  'is_BP_including_BPSCZ.pval' : 'is_BP_including_BPSCZ_pval',
  'is_BPNOS.pval' : 'is_BPNOS_pval',
  'is_BPSCZ.pval' : 'is_BPSCZ_pval',
  'is_BPPSY.pval' : 'is_BPPSY_pval',
  'is_BP_no_PSY.pval' : 'is_BP_no_PSY_pval',
  # Odds ratios
  'is_BP1.OR' : 'is_BP1_OR',
  'is_BP2.OR' : 'is_BP2_OR',
  'is_BP.OR' : 'is_BP_OR',
  'is_BP_including_BPSCZ.OR' : 'is_BP_including_BPSCZ_OR',
  'is_BPNOS.OR' : 'is_BPNOS_OR',
  'is_BPSCZ.OR' : 'is_BPSCZ_OR',
  'is_BPPSY.OR' : 'is_BPPSY_OR',
  'is_BP_no_PSY.OR' : 'is_BP_no_PSY_OR'
  }
)

ht = ht.annotate(**ht_fisher[ht.key])

# Next the Fisher's exact for MAC5 (not in gnomAD).
ht_fisher = hl.import_table("gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_fisher_tables/BP_gene_fisher_MAC5_gnom_non_psych_qq.tsv", impute=True).key_by('gene_symbol', 'consequence_category')
ht_fisher = ht_fisher.select(
  # p-values
  'gnom_non_psych.is_BP1.pval',
  'gnom_non_psych.is_BP2.pval',
  'gnom_non_psych.is_BP.pval',
  'gnom_non_psych.is_BP_including_BPSCZ.pval',
  'gnom_non_psych.is_BPNOS.pval',
  'gnom_non_psych.is_BPSCZ.pval',
  'gnom_non_psych.is_BPPSY.pval',
  'gnom_non_psych.is_BP_no_PSY.pval',
  # Odds ratios
  'gnom_non_psych.is_BP1.OR',
  'gnom_non_psych.is_BP2.OR',
  'gnom_non_psych.is_BP.OR',
  'gnom_non_psych.is_BP_including_BPSCZ.OR',
  'gnom_non_psych.is_BPNOS.OR',
  'gnom_non_psych.is_BPSCZ.OR',
  'gnom_non_psych.is_BPPSY.OR',
  'gnom_non_psych.is_BP_no_PSY.OR'
)

ht_fisher = ht_fisher.rename(
  {
  # p-values
  'gnom_non_psych.is_BP1.pval' : 'gnom_non_psych_is_BP1_pval',
  'gnom_non_psych.is_BP2.pval' : 'gnom_non_psych_is_BP2_pval',
  'gnom_non_psych.is_BP.pval' : 'gnom_non_psych_is_BP_pval',
  'gnom_non_psych.is_BP_including_BPSCZ.pval' : 'gnom_non_psych_is_BP_including_BPSCZ_pval',
  'gnom_non_psych.is_BPNOS.pval' : 'gnom_non_psych_is_BPNOS_pval',
  'gnom_non_psych.is_BPSCZ.pval' : 'gnom_non_psych_is_BPSCZ_pval',
  'gnom_non_psych.is_BPPSY.pval' : 'gnom_non_psych_is_BPPSY_pval',
  'gnom_non_psych.is_BP_no_PSY.pval' : 'gnom_non_psych_is_BP_no_PSY_pval',
  # Odds ratios
  'gnom_non_psych.is_BP1.OR' : 'gnom_non_psych_is_BP1_OR',
  'gnom_non_psych.is_BP2.OR' : 'gnom_non_psych_is_BP2_OR',
  'gnom_non_psych.is_BP.OR' : 'gnom_non_psych_is_BP_OR',
  'gnom_non_psych.is_BP_including_BPSCZ.OR' : 'gnom_non_psych_is_BP_including_BPSCZ_OR',
  'gnom_non_psych.is_BPNOS.OR' : 'gnom_non_psych_is_BPNOS_OR',
  'gnom_non_psych.is_BPSCZ.OR' : 'gnom_non_psych_is_BPSCZ_OR',
  'gnom_non_psych.is_BPPSY.OR' : 'gnom_non_psych_is_BPPSY_OR',
  'gnom_non_psych.is_BP_no_PSY.OR' : 'gnom_non_psych_is_BP_no_PSY_OR'
  }
)

ht = ht.annotate(**ht_fisher[ht.key])

# Next, the CMH tests
ht_CMH = hl.import_table("gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_CMH_tables/BP_gene_CMH_MAC5_qq.tsv",
  impute=True, missing=['', 'NA']).key_by('gene_symbol', 'consequence_category')
ht_CMH = ht_CMH.rename(
  {
  # p-values
  'is_BP1.pval' : 'is_BP1_CMH_pval',
  'is_BP2.pval' : 'is_BP2_CMH_pval',
  'is_BP.pval' : 'is_BP_CMH_pval',
  'is_BP_including_BPSCZ.pval' : 'is_BP_including_BPSCZ_CMH_pval',
  'is_BPNOS.pval' : 'is_BPNOS_CMH_pval',
  'is_BPSCZ.pval' : 'is_BPSCZ_CMH_pval',
  'is_BPPSY.pval' : 'is_BPPSY_CMH_pval',
  'is_BP_no_PSY.pval' : 'is_BP_no_PSY_CMH_pval',
  # Odds ratios
  'is_BP1.OR' : 'is_BP1_CMH_OR',
  'is_BP2.OR' : 'is_BP2_CMH_OR',
  'is_BP.OR' : 'is_BP_CMH_OR',
  'is_BP_including_BPSCZ.OR' : 'is_BP_including_BPSCZ_CMH_OR',
  'is_BPNOS.OR' : 'is_BPNOS_CMH_OR',
  'is_BPSCZ.OR' : 'is_BPSCZ_CMH_OR',
  'is_BPPSY.OR' : 'is_BPPSY_CMH_OR',
  'is_BP_no_PSY.OR' : 'is_BP_no_PSY_CMH_OR',
  }
)

ht = ht.annotate(**ht_CMH[ht.key])

ht_CMH = hl.import_table("gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_CMH_tables/BP_gene_CMH_MAC5_gnom_non_psych_qq.tsv",
  impute=True, missing=['', 'NA']).key_by('gene_symbol', 'consequence_category')
ht_CMH = ht_CMH.rename(
  {
  # p-values
  'gnom_non_psych.is_BP1.pval' : 'gnom_non_psych_is_BP1_CMH_pval',
  'gnom_non_psych.is_BP2.pval' : 'gnom_non_psych_is_BP2_CMH_pval',
  'gnom_non_psych.is_BP.pval' : 'gnom_non_psych_is_BP_CMH_pval',
  'gnom_non_psych.is_BP_including_BPSCZ.pval' : 'gnom_non_psych_is_BP_including_BPSCZ_CMH_pval',
  'gnom_non_psych.is_BPNOS.pval' : 'gnom_non_psych_is_BPNOS_CMH_pval',
  'gnom_non_psych.is_BPSCZ.pval' : 'gnom_non_psych_is_BPSCZ_CMH_pval',
  'gnom_non_psych.is_BPPSY.pval' : 'gnom_non_psych_is_BPPSY_CMH_pval',
  'gnom_non_psych.is_BP_no_PSY.pval' : 'gnom_non_psych_is_BP_no_PSY_CMH_pval',
  # Odds ratios
  'gnom_non_psych.is_BP1.OR' : 'gnom_non_psych_is_BP1_CMH_OR',
  'gnom_non_psych.is_BP2.OR' : 'gnom_non_psych_is_BP2_CMH_OR',
  'gnom_non_psych.is_BP.OR' : 'gnom_non_psych_is_BP_CMH_OR',
  'gnom_non_psych.is_BP_including_BPSCZ.OR' : 'gnom_non_psych_is_BP_including_BPSCZ_CMH_OR',
  'gnom_non_psych.is_BPNOS.OR' : 'gnom_non_psych_is_BPNOS_CMH_OR',
  'gnom_non_psych.is_BPSCZ.OR' : 'gnom_non_psych_is_BPSCZ_CMH_OR',
  'gnom_non_psych.is_BPPSY.OR' : 'gnom_non_psych_is_BPPSY_CMH_OR',
  'gnom_non_psych.is_BP_no_PSY.OR' : 'gnom_non_psych_is_BP_no_PSY_CMH_OR'
  }
)

ht = ht.annotate(**ht_CMH[ht.key])

ht = ht.annotate(explode_data = [
  hl.struct(analysis_group = 'Bipolar Disorder 1',
    case_count=ht.is_BP1.case_count, control_count=ht.is_BP1.control_count,
    n_cases=ht.n_is_BP1.cases, n_controls=ht.n_is_BP1.controls,
    # p-values
    fisher_pval= ht.is_BP1_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BP1_pval,
    CMH_pval= ht.is_BP1_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BP1_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BP1_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BP1_OR,
    CMH_OR= ht.is_BP1_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BP1_CMH_OR
    ),
  hl.struct(analysis_group = 'Bipolar Disorder 2',
    case_count=ht.is_BP2.case_count, control_count=ht.is_BP2.control_count,
    n_cases=ht.n_is_BP2.cases, n_controls=ht.n_is_BP2.controls,
    # p-values
    fisher_pval= ht.is_BP2_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BP2_pval,
    CMH_pval= ht.is_BP2_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BP2_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BP2_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BP2_OR,
    CMH_OR= ht.is_BP2_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BP2_CMH_OR
    ),
  # hl.struct(analysis_group = 'Bipolar Disorder NOS',
  #   case_count=ht.is_BPNOS.case_count, control_count=ht.is_BPNOS.control_count,
  #   n_cases=ht.n_is_BPNOS.cases, n_controls=ht.n_is_BPNOS.controls,
  #   # p-values
  #   fisher_pval= ht.is_BPNOS_pval,
  #   fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BPNOS_pval,
  #   CMH_pval= ht.is_BPNOS_CMH_pval,
  #   CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BPNOS_CMH_pval,
  #   # Odds ratios
  #   fisher_OR= ht.is_BPNOS_OR,
  #   fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BPNOS_OR,
  #   CMH_OR= ht.is_BPNOS_CMH_OR,
  #   CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BPNOS_CMH_OR,
  # ),
  # hl.struct(analysis_group = 'Schizoaffective',
  #   case_count=ht.is_BPSCZ.case_count, control_count=ht.is_BPSCZ.control_count,
  #   n_cases=ht.n_is_BPSCZ.cases, n_controls=ht.n_is_BPSCZ.controls,
  #   # p-values
  #   fisher_pval= ht.is_BPSCZ_pval,
  #   fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BPSCZ_pval,
  #   CMH_pval= ht.is_BPSCZ_CMH_pval,
  #   CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BPSCZ_CMH_pval,
  #   # Odds ratios
  #   fisher_OR= ht.is_BPSCZ_OR,
  #   fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BPSCZ_OR,
  #   CMH_OR= ht.is_BPSCZ_CMH_OR,
  #   CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BPSCZ_CMH_OR
  #   ),
  hl.struct(analysis_group = 'Bipolar Disorder',
    case_count=ht.is_BP.case_count, control_count=ht.is_BP.control_count,
    n_cases=ht.n_is_BP.cases, n_controls=ht.n_is_BP.controls,
    # p-values
    fisher_pval= ht.is_BP_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_pval,
    CMH_pval= ht.is_BP_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BP_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_OR,
    CMH_OR= ht.is_BP_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_CMH_OR
    ),
  hl.struct(analysis_group = 'Bipolar Disorder (including Schizoaffective)',
    case_count=ht.is_BP_including_BPSCZ.case_count, control_count=ht.is_BP_including_BPSCZ.control_count,
    n_cases=ht.n_is_BP_including_BPSCZ.cases, n_controls=ht.n_is_BP_including_BPSCZ.controls,
    # p-values
    fisher_pval= ht.is_BP_including_BPSCZ_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_including_BPSCZ_pval,
    CMH_pval= ht.is_BP_including_BPSCZ_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_including_BPSCZ_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BP_including_BPSCZ_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_including_BPSCZ_OR,
    CMH_OR= ht.is_BP_including_BPSCZ_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_including_BPSCZ_CMH_OR
    ),
  hl.struct(analysis_group = 'Bipolar Disorder with Psychosis',
    case_count=ht.is_BPPSY.case_count, control_count=ht.is_BPPSY.control_count,
    n_cases=ht.n_is_BPPSY.cases, n_controls=ht.n_is_BPPSY.controls,
    # p-values
    fisher_pval= ht.is_BPPSY_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BPPSY_pval,
    CMH_pval= ht.is_BPPSY_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BPPSY_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BPPSY_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BPPSY_OR,
    CMH_OR= ht.is_BPPSY_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BPPSY_CMH_OR
    ),
  hl.struct(analysis_group = 'Bipolar Disorder without Psychosis',
    case_count=ht.is_BP_no_PSY.case_count, control_count=ht.is_BP_no_PSY.control_count,
    n_cases=ht.n_is_BP_no_PSY.cases, n_controls=ht.n_is_BP_no_PSY.controls,
    # p-values
    fisher_pval= ht.is_BP_no_PSY_pval,
    fisher_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_no_PSY_pval,
    CMH_pval= ht.is_BP_no_PSY_CMH_pval,
    CMH_gnom_non_psych_pval= ht.gnom_non_psych_is_BP_no_PSY_CMH_pval,
    # Odds ratios
    fisher_OR= ht.is_BP_no_PSY_OR,
    fisher_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_no_PSY_OR,
    CMH_OR= ht.is_BP_no_PSY_CMH_OR,
    CMH_gnom_non_psych_OR= ht.gnom_non_psych_is_BP_no_PSY_CMH_OR
    )
  ]
)

ht = ht.select('explode_data')
ht = ht.explode('explode_data')
ht = ht.transmute(**ht.explode_data)

# Ensure that all floating point precision for taking logs is dealt with
# The following is when we used -log10p rather than p for displaying the results.
# ht = ht.transmute(
#   fisher_log_pval = hl.cond(ht.fisher_log_pval > 0, 0, ht.fisher_log_pval),
#   fisher_gnom_non_psych_log_pval = hl.cond(ht.fisher_gnom_non_psych_log_pval > 0, 0, ht.fisher_gnom_non_psych_log_pval),
#   CMH_log_pval = hl.cond(ht.CMH_log_pval > 0, 0, ht.CMH_log_pval),
#   CMH_gnom_non_psych_log_pval = hl.cond(ht.CMH_gnom_non_psych_log_pval > 0, 0, ht.CMH_gnom_non_psych_log_pval)
# )

# Grab the gene names from genenames.org
# curl -o hgnc.tsv 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=md_mim_id&col=md_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
# Move them to the cloud
# gsutil cp hgnc.tsv gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_names_for_browser/

# Read in and merge by gene name to get gene ID (which is ensembl ID).
ht_var_ann = hl.read_table(BROWSER_VARIANT_ANNOTATION_TABLE).key_by('gene_name')
ht = ht.annotate(gene_id = ht_var_ann[ht.gene_symbol].gene_id)

# Merge by ensembl ID to get gene description
ht_names = hl.import_table('gs://dalio_bipolar_w1_w2_hail_02/analysis/gene_names_for_browser/hgnc.tsv', impute=True).key_by('Ensembl ID(supplied by Ensembl)')
ht = ht.annotate(**ht_names[ht.gene_id])
ht = ht.rename({'Approved name' : 'gene_description'})

ht.write(BROWSER_GENE_RESULTS_TABLE, overwrite=True)
