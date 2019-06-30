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

for i in {0..499}
do
    if [ -e ${pathSequences}/fasta_parts/${suffix}_rm_part${i}.fa ]; then 
        if [ -e  ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out ]; then
	    export last=`tail -n 1 ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out | cut -f 1`
	    export tot=`grep -c ">" ${pathSequences}/${ref}/fasta_parts/${suffix}_rm_part${i}.fa`
	    export index=`grep ">" ${pathSequences}/${ref}/fasta_parts/${suffix}_rm_part${i}.fa | grep -n ${last} | cut -f 1 -d ':'`
	    
	    export ratio=$(($index * 100 / $tot))
	    
	    echo ${ref} ${tg} ${suffix} ${seq} "part" ${i} "index "${index}" out of "${tot}" "${ratio}"% done";
	else
	    echo "cannot find blast results"
	fi
    fi
done

#####################################################################
