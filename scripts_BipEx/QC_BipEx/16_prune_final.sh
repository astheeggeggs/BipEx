#!/bin/bash

reuse PLINK2

for i in `seq 1 22`;
do
	plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink_b38/final_qc.chr${i} \
	--indep 50 5 2 --out /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/16_final_qc.chr${i}
done  

plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink_b38/final_qc.chrX \
      --indep 50 5 2\
      --out /stanley/genetics/analysis/bipolar_dalio/plink_b38/17_final_qc.chrX

cd /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/
cat 16_final_qc.chr1.prune.in 16_final_qc.chr2.prune.in > 16_prune.final_qc.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat 16_prune.final_qc.keep.variant_list_tmp 16_final_qc.chr${i}.prune.in > 16_prune.final_qc.keep.variant_list
	mv 16_prune.final_qc.keep.variant_list 16_prune.final_qc.keep.variant_list_tmp
done

mv 16_prune.final_qc.keep.variant_list_tmp 16_prune.final_qc.keep.variant_list
gsutil cp 16_prune.final_qc.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/
