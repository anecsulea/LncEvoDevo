#!/bin/bash

export sp=$1
export annot=$2
export db=$3

#####################################################################

export path=LncEvoDevo
    
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathResults=${path}/results/coding_potential/${sp}
export pathScripts=${path}/scripts/coding_potential/blastx_proteins

#####################################################################

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathSequences=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathSequences=${pathStringTie}
    export suffix=FilteredTranscripts_StringTie_Ensembl${release}
fi

#####################################################################

if [ -e ${pathResults}/${suffix}_rm_vs_${db}.blastx.out ]; then
    echo "results file already here, not doing anything"
else    
    for i in {0..499}
    do
	if [ -e ${pathSequences}/fasta_parts/${suffix}_rm_part${i}.fa ]; then 
	    if [ -e ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out ]; then
		cat ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out >> ${pathResults}/${suffix}_rm_vs_${db}.blastx.out
	    else
		echo "weird! cannot find results for part "${i} ${sp} ${annot} ${db}
	    fi
	fi
    done
fi

#####################################################################
