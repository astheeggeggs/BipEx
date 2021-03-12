import hail as hl

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls_updated_phenotypes_rekeyed.mt'

def gene_burden_annotations(mt, bool_phenotype):

    mt = mt.group_rows_by(
        gene_symbol = mt.gene_symbol,
        vep_consequence = mt.vep_consequence,
        is_snp = mt.is_snp,
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))

    mt_results = mt.annotate_rows(
        counts = hl.agg.group_by(
            (mt.phenotype_boolean[bool_phenotype], mt.n > 0), hl.agg.count()
        )
    )

    ht = mt_results.rows()
    ht = ht.transmute(
        case_count = ht.counts.get((True, True), 0),
        case_no_count = ht.counts.get((True, False), 0),
        control_count = ht.counts.get((False, True), 0),
        control_no_count = ht.counts.get((False, False), 0)
    )

    return(ht)

def gene_burden_annotations_per_population(mt, bool_phenotype):

    mt = mt.group_rows_by(
        gene_symbol = mt.gene_symbol,
        vep_consequence = mt.vep_consequence,
        is_snp = mt.is_snp,
        ).aggregate(n = hl.agg.count_where(mt.GT.is_non_ref()))

    mt_results = mt.group_cols_by('forest_location').aggregate(
        counts = hl.agg.group_by(
            (mt.phenotype_boolean[bool_phenotype], mt.n > 0), hl.agg.count()
        )
    )

    ht = mt_results.entries()
    ht = ht.transmute(
        case_count = ht.counts.get((True, True), 0),
        case_no_count = ht.counts.get((True, False), 0),
        control_count = ht.counts.get((False, True), 0),
        control_no_count = ht.counts.get((False, False), 0)
    )

    return ht

def gene_burden_bool_phenotypes(mt, bool_phenotypes, annotation_grouping_function):
    
    init = True
    for bool_phenotype in bool_phenotypes:
        if init:
            ht_gene = annotation_grouping_function(mt, bool_phenotype)
            ht_gene = ht_gene.annotate(**{bool_phenotype:hl.struct()})
            ann_data = ht_gene[bool_phenotype].annotate(
                case_count = ht_gene.case_count,
                case_no_count = ht_gene.case_no_count,
                control_count = ht_gene.control_count,
                control_no_count = ht_gene.control_no_count
            )
            ht_gene = ht_gene.transmute(**{bool_phenotype: ann_data})
        else:
            ht_gene_tmp = annotation_grouping_function(mt, bool_phenotype)
            ht_gene = ht_gene.annotate(**{bool_phenotype:hl.struct()})
            ann_data = ht_gene[bool_phenotype].annotate(
                case_count = ht_gene_tmp[ht_gene.key].case_count,
                case_no_count = ht_gene_tmp[ht_gene.key].case_no_count,
                control_count = ht_gene_tmp[ht_gene.key].control_count,
                control_no_count = ht_gene_tmp[ht_gene.key].control_no_count
            )
            ht_gene = ht_gene.annotate(**{bool_phenotype: ann_data})
        init = False

    return(ht_gene)

def gene_burden_add_case_control_counts(ht_gene, mt, bool_phenotypes):

    for bool_phenotype in bool_phenotypes:
        ht_gene = ht_gene.annotate_globals(**{'n_' + bool_phenotype:hl.struct()})
        ann_data = ht_gene['n_' + bool_phenotype].annotate(
            cases = mt.aggregate_cols(hl.agg.sum(mt['phenotype_boolean'][bool_phenotype])),
            controls = mt.aggregate_cols(hl.agg.sum(~mt['phenotype_boolean'][bool_phenotype]))
        )
        ht_gene = ht_gene.annotate_globals(**{'n_' + bool_phenotype:ann_data})

    return(ht_gene)

def gene_burden_bool_phenotypes_all_classes(mt, bool_phenotypes, annotation_grouping_function):

    ht_gene = gene_burden_bool_phenotypes(mt, bool_phenotypes, annotation_grouping_function)
    ht_gene = ht_gene.annotate(
        burden = ht_gene[ht_gene.key]
        )
    mt = mt.filter_rows(~mt.annotation.inGnomAD_nonpsych)
    ht_gene = ht_gene.annotate(
        burden_gnom_non_psych = gene_burden_bool_phenotypes(
            mt, bool_phenotypes, annotation_grouping_function)[ht_gene.key]
        )
    
    # Add global annotations, restrict to the variables that we require, and write to .tsv and ht.
    ht_gene = gene_burden_add_case_control_counts(ht_gene, mt, bool_phenotypes)
    
    return(ht_gene)

# We require a function that filters to the gene and then performs all the counts. 
mt = hl.read_matrix_table(QC_HARDCALLS_SHUFFLE_AVOIDANCE_MT)

# Gene to look at
gene_name = "AKAP11"

mt = mt.filter_rows(mt.gene_symbol == gene_name)
mt = mt.annotate_rows(
  vep_consequence = mt.annotation.vep.worst_csq_for_variant_canonical.most_severe_consequence
  )
mt = mt.annotate_rows(is_snp = hl.is_snp(mt.alleles[0], mt.alleles[1]))
mt_MAC5 = mt.filter_rows(mt.is_MAC5)

bool_phenotypes = ['is_BP1', 'is_BP2', 'is_BP', 'is_BPPSY', 'is_BP_no_PSY']
mt_BP = mt_MAC5.filter_cols(hl.is_defined(mt_MAC5['phenotype_boolean']['is_BP_including_BPSCZ']))

ht_gene_BP = gene_burden_bool_phenotypes_all_classes(mt_BP, bool_phenotypes, gene_burden_annotations)
ht_gene_BP.show()
