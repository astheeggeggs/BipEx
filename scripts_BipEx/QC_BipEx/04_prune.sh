#!/bin/bash

mkdir /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned

reuse PLINK2

for i in `seq 1 22`;
do
	plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink_b38/filterGT.chr${i} \
	--indep 50 5 2 --out /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/04_chr${i}
done 

plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink_b38/filterGT.chrX \
      --indep 50 5 2 --out /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/04_chrX

gsutil cp /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/04_chrX.prune.in gs://dalio_bipolar_w1_w2_hail_02/data/variants/

cd /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned
cat 04_chr1.prune.in 04_chr2.prune.in > 04_prune.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat 04_prune.keep.variant_list_tmp 04_chr${i}.prune.in > 04_prune.keep.variant_list
	mv 04_prune.keep.variant_list 04_prune.keep.variant_list_tmp
done

mv 04_prune.keep.variant_list_tmp 04_prune.keep.variant_list
gsutil cp 04_prune.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/
