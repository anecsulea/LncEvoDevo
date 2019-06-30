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

###################################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_splice

for sample in `ls ${pathResults} | grep -v txt `
do
    if [ -e ${pathResults}/${sample}/transcripts_splice_junctions.txt ]; then
	echo "already done"
    else
	echo "perl ${pathScripts}/check.splice.junctions.pl --pathAnnotGTF=${pathResults}/${sample}/${filename}.gtf --pathCorrectJunctions=${pathAlignments}/${sample}/junctions.txt --pathWrongJunctions=${pathAlignments}/${sample}/junctions_wrongstrand.txt  --pathOutput=${pathResults}/${sample}/transcripts_splice_junctions.txt"  >>  ${pathScripts}/bsub_script_splice
    fi
done

qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_splice_${sp}.txt -e ${pathScripts}/std_error_splice_${sp}.txt  ${pathScripts}/bsub_script_splice


###################################################################################
