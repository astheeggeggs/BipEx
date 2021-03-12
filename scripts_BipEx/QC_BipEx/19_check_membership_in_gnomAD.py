import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Want to check whether any of the samples are in gnomAD
# Check distribution of singletons not in gnomAD - should be non-zero!

QC_HARDCALLS_MT = 'gs://dalio_bipolar_w1_w2_hail_02/data/mt/17_european.strict.hardcalls.mt'
GNOMAD_TSV = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/19_gnomAD_check.tsv'

mt = hl.read_matrix_table(QC_HARDCALLS_MT)

mt = mt.filter_cols((mt.phenotype.PHENOTYPE_COARSE == "Bipolar Disorder") | (mt.phenotype.PHENOTYPE_COARSE == "Control"))
mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
mt = mt.filter_rows(mt.is_singleton)
mt = mt.annotate_cols(singleton_count = hl.agg.count_where(mt.GT.is_non_ref()))

mt = mt.filter_rows((~mt.annotation.inGnomAD_nonpsych) & (mt.is_singleton))
mt = mt.annotate_cols(not_inGnomAD_count = hl.agg.count_where(mt.GT.is_non_ref()))

scatter = hl.plot.scatter(mt.not_inGnomAD_count, mt.singleton_count)
show(scatter)

mt.cols().select('phenotype', 'singleton_count', 'not_inGnomAD_count').flatten().export(GNOMAD_TSV)
