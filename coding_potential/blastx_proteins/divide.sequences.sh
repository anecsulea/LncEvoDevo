#!/bin/bash

export sp=$1
export annot=$2
export nbparts=$3

#####################################################################

export path=LncEvoDevo
	
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathScripts=${path}/scripts/coding_potential/blastx_proteins

export release=94

###################################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathSequences=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathSequences=${pathStringTie}
    export suffix=FilteredTranscripts_StringTie_Ensembl${release}
fi

###################################################################################

if [ -e ${pathSequences}/fasta_parts ]; then
    echo "path output already there"
else
   mkdir ${pathSequences}/fasta_parts
fi

###################################################################################

perl ${pathScripts}/divide.sequences.pl --pathFastaInput=${pathSequences}/${suffix}_rm.fa --nbParts=${nbparts} --prefixOutput=${pathSequences}/fasta_parts/${suffix}_rm

###################################################################################

