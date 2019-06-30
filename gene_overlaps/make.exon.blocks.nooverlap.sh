#!/bin/bash

export species=$1

#################################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathGeneOverlaps=${path}/results/gene_overlaps/${species}
export pathScripts=${path}/scripts/gene_overlaps

export release=94

export prefix=FilteredTranscripts_StringTie_Ensembl${release}

#################################################################################

perl make.exon.blocks.nooverlap.pl --pathExonCoords=${pathGeneOverlaps}/ExonCoordinates_ExcludingOverlapOtherGenes_${prefix}.txt --collapseDistance=0 --pathOutput=${pathGeneOverlaps}/ExonBlocks_ExcludingOverlapOtherGenes_${prefix}.txt

#################################################################################

    
