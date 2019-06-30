#!/bin/bash

######################################################

export path=LncEvoDevo
export pathStringTie=${path}/results/stringtie_assembly
export pathDatasets=${path}/supplementary_datasets

######################################################

export release=94

######################################################

for sp in Mouse Rat Chicken
do
    cp ${pathStringTie}/${sp}/combined/FilteredTranscripts_StringTie_Ensembl${release}.gtf ${pathDatasets}/SupplementaryDataset1/FullAnnotation_${sp}.gtf
done

######################################################

