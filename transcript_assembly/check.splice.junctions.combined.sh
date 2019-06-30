#!/bin/bash

export sp=$1
export type=$2

###################################################################################

export path=LncEvoDevo

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

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_splice

export pathCorrectJunctions=""
export pathWrongJunctions=""

for sample in `ls ${pathAlignments} | grep -v txt `
do
    export pathCorrectJunctions=${pathAlignments}/${sample}/junctions.txt,${pathCorrectJunctions}
    export pathWrongJunctions=${pathAlignments}/${sample}/junctions_wrongstrand.txt,${pathWrongJunctions}
done

echo "perl ${pathScripts}/check.splice.junctions.pl --pathAnnotGTF=${pathResults}/combined/FilteredTranscripts_StringTie_Ensembl${release}.gtf --pathCorrectJunctions=${pathCorrectJunctions} --pathWrongJunctions=${pathWrongJunctions} --pathOutput=${pathResults}/combined/FilteredTranscripts_StringTie_Ensembl${release}_SpliceJunctionsStats.txt"  >>  ${pathScripts}/bsub_script_splice

qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_splice_${sp}.txt -e ${pathScripts}/std_error_splice_${sp}.txt  ${pathScripts}/bsub_script_splice

###################################################################################
