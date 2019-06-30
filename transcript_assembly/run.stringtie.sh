#!/bin/bash

export species=$1

#############################################################################

export path=LncEvoDevo
    
export pathDocs=${path}/docs
export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathAlignments=${path}/results/hisat/${species}
export pathResults=${path}/results/stringtie_assembly/${species}
export pathScripts=${path}/scripts/transcript_assembly

export release=94

#############################################################################

for sample in `ls ${pathAlignments} | grep -v txt`
do
    export strand="--rf"
    
    if [ ${species} = "Chicken" ]; then
	export library=`grep ^${sample}$'\t' ${pathDocs}/LibraryType_Chicken.txt | cut -f 2`
	
	if [ ${library} = "fr-unstranded" ]; then
	    export strand=""
	fi
	
	if [ ${library} = "fr-firststrand" ]; then
	    export strand="--rf"
	fi
	
	if [ ${library} = "fr-secondstrand" ]; then
	    export strand="--fr"
	fi
    fi

    if [ -e ${pathResults}/${sample} ]; then
	echo "output already there"
    else
	mkdir ${pathResults}/${sample}
    fi

    if [ -e ${pathResults}/${sample}/assembled_transcripts.gtf ]&&[ -e ${pathResults}/${sample}/gene_abundance.txt ]; then
	echo "already done"
    else
	if [ -e ${pathAlignments}/${sample}/accepted_hits.bam ]; then
	  
	    echo "#!/bin/bash" > ${pathScripts}/bsub_script_stringtie
	   
    	    echo "stringtie ${pathAlignments}/${sample}/accepted_hits.bam -G ${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf -m 150 -a 8 -f 0.05 -p 8 -o ${pathResults}/${sample}/assembled_transcripts.gtf -A ${pathResults}/${sample}/gene_abundance.txt ${strand}">> ${pathScripts}/bsub_script_stringtie
	    
	    qsub -P P_biometr -q mc_huge -l s_rss=8G,sps=1 -pe multicores 8 -o ${pathScripts}/std_output_stringtie_${sp}.txt -e ${pathScripts}/std_error_stringtie_${sp}.txt ${pathScripts}/bsub_script_stringtie
	    
	fi
    fi
    
done

#############################################################################
