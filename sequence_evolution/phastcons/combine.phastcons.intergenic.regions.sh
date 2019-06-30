#!/bin/bash

export species=$1
export annot=$2
export phast=$3

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathResults=${path}/results/sequence_evolution/phastcons/${species}/${phast}
export pathPhastCons=${path}/data/phastcons/${species}

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}
    export suffixExons=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}
    export suffixExons=FilteredTranscripts_StringTie_Ensembl${release}
fi


export minsize=1kb
export mindist=5kb

#####################################################################

cp ${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_chr1_${suffixExons}.txt ${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_${suffixExons}.txt 

for chr in {2..19} X Y
do
    sed '1d'  ${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_chr${chr}_${suffixExons}.txt >> ${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_${suffixExons}.txt 
done

#####################################################################
