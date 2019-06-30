#!/bin/bash

export species=$1

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathEnsembl=${path}/data/genome_indexes/${species}
export pathUCSC=${path}/data/UCSC_sequences/${species}
export pathScripts=${path}/scripts/exon_projections

export release=87

####################################################################################

if [ ${species} = "Mouse" ]; then
    export suffix=mm10
fi

if [ ${species} = "Rat" ]; then
    export suffix=rn6
fi

if [ ${species} = "Chicken" ]; then
    export suffix=galGal4
fi


if [ ${species} = "Human" ]; then
    export release=89
    export suffix=hg38
fi

####################################################################################

if [ -e ${pathEnsembl}/genome_ensembl${release}.fa ]; then
    perl ${pathScripts}/chromosome.correspondence.pl --pathEnsembl=${pathEnsembl}/genome_ensembl${release}.fa --pathUCSC=${pathUCSC}/${suffix}.fa.gz --pathOutput=${pathUCSC}/chromosomes_Ensembl_UCSC.txt
else
    if [ -e ${pathEnsembl}/genome_ensembl${release}.fa.gz ]; then
	perl ${pathScripts}/chromosome.correspondence.pl --pathEnsembl=${pathEnsembl}/genome_ensembl${release}.fa.gz --pathUCSC=${pathUCSC}/${suffix}.fa.gz --pathOutput=${pathUCSC}/chromosomes_Ensembl_UCSC.txt
    fi
fi

####################################################################################
