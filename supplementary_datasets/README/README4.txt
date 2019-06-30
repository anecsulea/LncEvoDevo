#####################################################################

This supplementary dataset contains the results of the differential expression analyses, testing for temporal variations in protein-coding gene and lncRNA expression. The analyses were performed separately for each organ and the results were grouped into a single file for all organs. The analyses looking for a global age effect  were done only for mouse and rat, as for chicken we only sampled two developmental stages.

DifferentialExpression_GlobalAgeEffect_AllReads_${species}.txt: differential expression analysis, testing for a global effect of the developmental stage. All available reads were used for this analysis.

DifferentialExpression_GlobalAgeEffect_ResampledReads_${species}.txt: differential expression analysis, testing for a global effect of the developmental stage, after downsampling protein-coding genes read counts to match those of lncRNAs.

DifferentialExpression_ConsecutiveStages_${species}.txt: diffential expression analysis, testing for a difference between two consecutive developmental stages.

KmeansClusters_MouseRatOrtho_DiffExp_${organ}_ProteinCoding.txt: results of K-means clustering analysis for protein-coding genes that are differentially expressed among developmental stages, for each organ, in both mouse and rat. These results are shown in Supplementary Figure 13.

GOEnrichment_BiologicalProcess_KmeansClusters_MouseRatOrtho_DiffExp_${organ}_Cluster${i}_ProteinCoding.txt: results of a gene ontology enrichment analysis, for each K-means cluster defined above.


KmeansClusters_MouseRatOrtho_DiffExp_${organ}_LncRNAs.txt: results of K-means clustering analysis for lncRNAs that are differentially expressed among developmental stages, for each organ, in both mouse and rat. These results are shown in Figure 8.


#####################################################################

GOEnrichment_DifferentialExpression_ConsecutiveStages subfolder:

BiologicalProcess_${direction}regulatedGenes_${species}_${organ}_${stage1}_${stage2}.txt: gene ontology enrichment analysis for protein-coding genes that are up- or down-regulated in a given species and organ, between two consecutive developmental stages.
