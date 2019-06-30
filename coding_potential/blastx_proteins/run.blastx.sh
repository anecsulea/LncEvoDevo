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

if [ ${db} = "SwissProt" ]; then
    export pathBlastDB=${path}/data/SwissProt/uniprot_sprot_filtered
fi

if [ ${db} = "Pfam" ]; then
    export pathBlastDB=${path}/data/Pfam/Pfam-A
fi

#####################################################################

if [ -e ${pathResults}/blastx_parts ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}/blastx_parts
fi

#####################################################################

for i in {0..499}
do
    if [ -e ${pathSequences}/fasta_parts/${suffix}_rm_part${i}.fa ]; then 
	if [ -e ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out ]; then
	    echo "already done"
	else
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_blast
	   
	    echo "blastx -num_threads 4 -query ${pathSequences}/fasta_parts/${suffix}_rm_part${i}.fa -db ${pathBlastDB} -out ${pathResults}/blastx_parts/${suffix}_rm_vs_${db}_part${i}.blastx.out -outfmt 6 -max_target_seqs 20 -seg yes -strand plus" >> ${pathScripts}/bsub_script_blast
	    
	    qsub -q mc_huge -l s_rss=4G,sps=1 -pe multicores 4 -o ${pathScripts}/std_output_blast_${sp}_${annot}_${db}.txt -e ${pathScripts}/std_error_blast_${sp}_${annot}_${db}.txt ${pathScripts}/bsub_script_blast
	   
	fi
    fi
done

#####################################################################
