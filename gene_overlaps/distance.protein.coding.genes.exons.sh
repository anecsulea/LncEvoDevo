#!/bin/bash

export sp=$1

#####################################################################

export path=LncEvoDevo
   
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#####################################################################

perl ${pathScripts}/distance.protein.coding.genes.exons.pl --pathGTF1=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathGTF2=${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=NA  --forbiddenTranscripts=ENSMUST00000127091,ENSMUST00000127664 --pathOutput=${pathResults}/Distance_ProteinCodingGeneExons_SenseStrand_AllTranscripts_Ensembl${release}.txt 

#####################################################################

perl ${pathScripts}/distance.protein.coding.genes.exons.pl --pathGTF1=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathGTF2=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=NA  --forbiddenTranscripts=ENSMUST00000127091,ENSMUST00000127664 --pathOutput=${pathResults}/Distance_ProteinCodingGeneExons_SenseStrand_FilteredTranscripts_StringTie_Ensembl${release}.txt 

#####################################################################
