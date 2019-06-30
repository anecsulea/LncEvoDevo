#!/bin/bash

export sp=$1
export annot=$2

###################################################################################

export path=LncEvoDevo
   
export pathAlignments=${path}/results/hisat/${sp}
export pathGenomeBrowser=${path}/results/genome_browser/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathAnnot=${path}/data/ensembl_annotations/${sp}
export pathResults=${path}/results/antisense_transcription/${sp}
export pathScripts=${path}/scripts/antisense_transcription

export release=94

###############################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathAnnot}/FilteredTranscripts_Ensembl${release}.gtf
    export suffixOutput=FilteredTranscripts_Ensembl${release}.txt
fi

if [ ${annot} = "StringTie" ]; then
    export pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf
    export suffixOutput=FilteredTranscripts_StringTie_Ensembl${release}.txt
fi

###################################################################################

for sample in `ls ${pathAlignments} | grep -v txt ` 
do
    if [ -e ${pathResults}/${sample} ]; then
	echo "path output already there"
    else
	mkdir ${pathResults}/${sample} 
    fi

    if [ -e ${pathResults}/${sample}/CoverageExons_${suffixOutput} ]&&[ -e ${pathResults}/${sample}/CoverageTranscripts_${suffixOutput} ]; then
	echo "already done"
    else
	if [ -e ${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz ]&&[ -e ${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz ]; then
	    
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage

	    echo "perl ${pathScripts}/compute.transcript.coverage.pl --pathAnnotGTF=${pathGTF} --pathCoverageForward=${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz --pathCoverageReverse=${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz --pathOutputExons=${pathResults}/${sample}/CoverageExons_${suffixOutput} --pathOutputTranscripts=${pathResults}/${sample}/CoverageTranscripts_${suffixOutput}"  >>  ${pathScripts}/bsub_script_coverage
	   
	    qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_coverage_${sample}_${sp}_${annot}.txt -e ${pathScripts}/std_error_coverage_${sample}_${sp}_${annot}.txt ${pathScripts}/bsub_script_coverage
	   

	fi
    fi	    
done

###################################################################################
