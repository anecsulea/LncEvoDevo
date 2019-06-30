########################################################################

This supplementary dataset contains the results of the mutual information network analyis performed with ARACNE-AP.

FullNetwork_${species}.txt.gz: these files contain all predicted interactions between lncRNAs and other genes (including both lncRNAs and protein-coding genes), with orthologues in mouse and rat.

FilteredNetwork_FDR0.001_${species}.txt.gz: these files contain filtered predicted interactions, after applying an FDR threshold of 0.001.

SharedInteractions_FDR0.001.txt: this file contains interactions that are found in both mouse and rat.

GOEnrichment_BiologicalProcess_lncRNAs_Min50ConnectedGenes.txt.gz: this file contains the results of GO enrichment analysis, performed for lncRNAs that are connected with at least 50 protein-coding genes in both mouse and rat. The enrichment analysis contrasts the shared connecting protein-coding genes, with the entire set of 1-to-1 orthologous protein-coding genes used for the mutual information analysis.

########################################################################
