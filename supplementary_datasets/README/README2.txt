#####################################################################

This supplementary dataset contains gene expression levels, measured with different approaches.

- KallistoRawTPM_${species}.txt: Kallisto TPM values  for each gene.
- KallistoNormalizedTPM_${species}.txt: Kallisto TPM values, after normalization for between-samples comparison with the median scaling procedure described in Brawand et al., Nature, 2011.
- KallistoEstimatedCounts_${sp}.txt: Kallisto estimated read counts for each gene. 
- UniqueReadCounts_${sp}.txt: unique read counts for each gene, computed with featureCounts in the Subread R package.
- UniqueReadCounts_Downsampled${nbreads}_${species}.txt: unique read count for each gene, after downsampling the same number of uniquely mapped reads for each organ/developmental stage combination.

#####################################################################

