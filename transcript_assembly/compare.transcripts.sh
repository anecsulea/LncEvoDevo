#!/bin/bash

export sp=$1
export type=$2

###################################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathAlignments=${path}/results/hisat/${sp}
export pathGenomeBrowser=${path}/results/genome_browser/${sp}
export pathResults=${path}/results/${type}_assembly/${sp}
export pathScripts=${path}/scripts/transcript_assembly


if [ ${type} = "cufflinks" ]; then
    export filename="transcripts"
else
    if [ ${type} = "stringtie" ]; then
	export filename="assembled_transcripts"
    else
	echo "unknown type"
	exit
    fi
fi

export release=94

###################################################################################
 
echo "#!/bin/bash" >  ${pathScripts}/bsub_script_comparetx

for sample in `ls ${pathResults} | grep -v txt | grep -v combined`
do   
    echo "perl ${pathScripts}/compare.transcripts.pl --pathAssembledGTF=${pathResults}/${sample}/${filename}.gtf --pathKnownGTF=${pathEnsembl}/FilteredTranscripts_Ensembl${release}.gtf  --pathOutput=${pathResults}/${sample}/comparison_Ensembl${release}.txt"  >>  ${pathScripts}/bsub_script_comparetx
done

qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_comparetx_${sp}.txt -e ${pathScripts}/std_error_comparetx_${sp}.txt ${pathScripts}/bsub_script_comparetx
 
 
###################################################################################
