#!/bin/bash

export sp=$1
export annot=$2

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/gene_overlaps/${sp}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export filename=FilteredTranscripts_Ensembl${release}
fi 


if [ ${annot} = "StringTie" ]; then
    export pathAnnot=${pathStringTie}
    export filename=FilteredTranscripts_StringTie_Ensembl${release}
fi 

#####################################################################

perl ${pathScripts}/extract.intron.blocks.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${filename}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=NA  --biotypes=protein_coding --pathOutputIntronBlocks=${pathAnnot}/IntronBlocks_ProteinCodingGenes_${filename}.txt --pathOutputSelectedExonBlocks=${pathAnnot}/ExonBlocks_ProteinCodingGenes_${filename}.txt  


perl ${pathScripts}/extract.intron.blocks.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${filename}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=NA --biotypes=all --pathOutputIntronBlocks=${pathAnnot}/IntronBlocks_AllGenes_${filename}.txt 

#####################################################################
