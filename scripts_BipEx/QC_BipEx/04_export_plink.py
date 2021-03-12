import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
PLINK_FILES = 'gs://dalio_bipolar_w1_w2_hail_02/data/plink/filterGT'

HIGH_LD_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/b38_high_ld.bed'
INITIAL_SAMPLES = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

high_LD_intervals = hl.import_locus_intervals(HIGH_LD_INTERVALS, reference_genome='GRCh38')

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.annotate_rows(in_high_LD = hl.is_defined(high_LD_intervals[mt.locus]))

mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]) & (~mt.in_high_LD))
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.01) & (mt.qc.AF [0]< 0.99) & ((mt.qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())).persist()

mt.count()

# def rename_samples(mt, mapping):
#     return mt.key_cols_by(s = hl.literal(mapping).get(mt.s, default=mt.s))

# mt = rename_samples(mt, {'431-BG00852 D':'431-BG00852_D'})

for x in range(1,23):

	mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval(hl.eval('chr' + hl.str(x)), reference_genome='GRCh38')])
	n_chr = mt_chr.count_rows()

	print('\nn variants in chr')
	print(x)
	print(n_chr)

	hl.export_plink(mt_chr, PLINK_FILES + '.chr' + str(x))

mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX', reference_genome='GRCh38')])
n_chr = mt_chr.count_rows()

print('\nn variants in chrX')
print(n_chr)

hl.export_plink(mt_chr, PLINK_FILES + '.chr' + 'X')
