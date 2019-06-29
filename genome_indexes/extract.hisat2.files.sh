#!/bin/bash 

export sp=$1
export release=$2

#####################################################################

export path=LncEvoDevo
export pathHisat2=${path}/data/genome_indexes/${sp}
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathScripts=${path}/scripts/genome_indexes

## hisat 2.0.5

######################################################################

hisat2_extract_splice_sites.py ${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf > ${pathEnsembl}/SpliceSites_Ensembl${release}.txt

hisat2_extract_exons.py ${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf > ${pathEnsembl}/Exons_Ensembl${release}.txt

######################################################################
