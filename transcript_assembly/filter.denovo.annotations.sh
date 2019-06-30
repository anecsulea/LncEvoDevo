#!/bin/bash

export species=$1
export type=$2

###########################################################################

export path=LncEvoDevo

###########################################################################

export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathResults=${path}/results/${type}_assembly/${species}
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

###########################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_filter


for sample in `ls ${pathResults} | grep -v txt | grep -v combined`
do
        
    echo "perl ${pathScripts}/filter.denovo.annotations.pl --pathAssembledGTF=${pathResults}/${sample}/${filename}.gtf --pathEnsemblSynonyms=${pathResults}/${sample}/comparison_Ensembl${release}.txt --pathReadthroughTranscripts=${pathResults}/${sample}/readthrough_transcripts.txt --minRatioSenseAntisense=0.01 --pathCoverage=${pathResults}/${sample}/coverage_transcripts.txt --pathJunctions=${pathResults}/${sample}/transcripts_splice_junctions.txt --pathOutputGTF=${pathResults}/${sample}/filtered_transcripts.gtf --pathOutputDiscardedTranscripts=${pathResults}/${sample}/discarded_transcripts.txt" >>  ${pathScripts}/bsub_script_filter
    
done

qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_filter_${species}.txt -e ${pathScripts}/std_error_filter_${species}.txt ${pathScripts}/bsub_script_filter

 
###########################################################################

