#!/bin/bash

##########################################################

export path=LncEvoDevo
export pathEnsemblOrtho=${path}/data/ensembl_ortho
export pathOurOrtho=${path}/results/lncRNA_dataset
export pathDatasets=${path}/supplementary_datasets

export release=94

##########################################################

cp ${pathEnsemblOrtho}/GeneFamilies_1to1_MouseRat_ProteinCoding_Ensembl${release}.txt ${pathDatasets}/SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Rat.txt
cp ${pathEnsemblOrtho}/GeneFamilies_1to1_MouseChicken_ProteinCoding_Ensembl${release}.txt ${pathDatasets}/SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Chicken.txt
cp ${pathEnsemblOrtho}/GeneFamilies_1to1_RatChicken_ProteinCoding_Ensembl${release}.txt ${pathDatasets}/SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Rat_Chicken.txt
cp ${pathEnsemblOrtho}/GeneFamilies_1to1_MouseRatChicken_ProteinCoding_Ensembl${release}.txt ${pathDatasets}/SupplementaryDataset5/EnsemblOrtho_ProteinCodingGenes_1to1_Mouse_Rat_Chicken.txt

cp ${pathOurOrtho}/GeneFamilies_FilteredTranscripts_StringTie_Ensembl${release}.txt ${pathDatasets}/SupplementaryDataset5/PredictedOrthoFamilies_Mouse_Rat_Chicken.txt

##########################################################
