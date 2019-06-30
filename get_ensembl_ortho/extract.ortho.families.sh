#!/bin/bash

##########################################################################

export path=LncEvoDevo

####################################################################################

export pathOrtho=${path}/data/ensembl_ortho
export pathAnnot=${path}/data/ensembl_annotations
export pathScripts=${path}/scripts/get_ensembl_ortho

export release=94

##########################################################################

# perl ${pathScripts}/extract.all1to1.ortho.families.pl --pathHomology=${pathOrtho}/homology_members_one2one_ensembl${release}.txt --pathEnsemblIDs=${pathOrtho}/EnsemblIDs.txt --speciesList=Mouse,Rat --prefixGeneInfo=${pathAnnot}/ --suffixGeneInfo=GeneInfo_Ensembl${release}.txt --pathOutputFamilies=${pathOrtho}/GeneFamilies_1to1_MouseRat_ProteinCoding_Ensembl${release}.txt

# perl ${pathScripts}/extract.all1to1.ortho.families.pl --pathHomology=${pathOrtho}/homology_members_one2one_ensembl${release}.txt --pathEnsemblIDs=${pathOrtho}/EnsemblIDs.txt --speciesList=Mouse,Rat,Chicken --prefixGeneInfo=${pathAnnot}/ --suffixGeneInfo=GeneInfo_Ensembl${release}.txt --pathOutputFamilies=${pathOrtho}/GeneFamilies_1to1_MouseRatChicken_ProteinCoding_Ensembl${release}.txt

##########################################################################

perl ${pathScripts}/extract.all1to1.ortho.families.pl --pathHomology=${pathOrtho}/homology_members_one2one_ensembl${release}.txt --pathEnsemblIDs=${pathOrtho}/EnsemblIDs.txt --speciesList=Mouse,Chicken --prefixGeneInfo=${pathAnnot}/ --suffixGeneInfo=GeneInfo_Ensembl${release}.txt --pathOutputFamilies=${pathOrtho}/GeneFamilies_1to1_MouseChicken_ProteinCoding_Ensembl${release}.txt

perl ${pathScripts}/extract.all1to1.ortho.families.pl --pathHomology=${pathOrtho}/homology_members_one2one_ensembl${release}.txt --pathEnsemblIDs=${pathOrtho}/EnsemblIDs.txt --speciesList=Rat,Chicken --prefixGeneInfo=${pathAnnot}/ --suffixGeneInfo=GeneInfo_Ensembl${release}.txt --pathOutputFamilies=${pathOrtho}/GeneFamilies_1to1_RatChicken_ProteinCoding_Ensembl${release}.txt

##########################################################################
