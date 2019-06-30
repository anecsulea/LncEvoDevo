#!/bin/bash

export species=$1
export release=$2

#################################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${species}

#################################################################################

perl make.exon.blocks.ensembl.pl --pathExonCoords=${pathEnsembl}/ExonCoords_Ensembl${release}.txt --pathExonAssignment=${pathEnsembl}/ExonsTranscripts_Ensembl${release}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathTranscriptInfo=${pathEnsembl}/TranscriptInfo_Ensembl${release}.txt --collapseDistance=10 --filter="yes" --pathReadThrough=${pathEnsembl}/ReadthroughTranscripts_Ensembl${release}.txt  --pathOutputExonBlocks=${pathEnsembl}/ExonBlocks_FilteredTranscripts_Ensembl${release}.txt

perl make.exon.blocks.ensembl.pl --pathExonCoords=${pathEnsembl}/ExonCoords_Ensembl${release}.txt --pathExonAssignment=${pathEnsembl}/ExonsTranscripts_Ensembl${release}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathTranscriptInfo=${pathEnsembl}/TranscriptInfo_Ensembl${release}.txt --collapseDistance=10 --filter="no" --pathReadThrough=NA --pathOutputExonBlocks=${pathEnsembl}/ExonBlocks_AllTranscripts_Ensembl${release}.txt

#################################################################################

    
