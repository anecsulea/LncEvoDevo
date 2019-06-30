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

for sample in `ls ${pathResults} | grep -v txt ` 
do
    if [ -e ${pathResults}/${sample}/coverage_exons.txt ]&&[ -e ${pathResults}/${sample}/coverage_transcripts.txt ]; then
	echo "already done"
    else
	if [ -e ${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz ]&&[ -e ${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz ]; then
	    
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage

	    echo "perl ${pathScripts}/compute.transcript.coverage.pl --pathAnnotGTF=${pathResults}/${sample}/${filename}.gtf --pathCoverageForward=${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz --pathCoverageReverse=${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz --pathOutputExons=${pathResults}/${sample}/coverage_exons.txt  --pathOutputTranscripts=${pathResults}/${sample}/coverage_transcripts.txt"  >>  ${pathScripts}/bsub_script_coverage
	    
	    qsub -q huge -l s_rss=6G,sps=1 -o ${pathScripts}/std_output_coverage_${sample}_${sp}.txt -e ${pathScripts}/std_error_coverage_${sample}_${sp}.txt ${pathScripts}/bsub_script_coverage
	    

	fi
    fi	    
done

###################################################################################
