#!/bin/bash

export sp=$1

######################################################################################################

export path=LncEvoDevo

export pathHisat=${path}/results/hisat/${sp}
export pathGenomeBrowser=${path}/results/genome_browser/${sp}
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathResults=${path}/results/RNA_degradation/${sp}
export pathScripts=${path}/scripts/RNA_degradation

export release=94

######################################################################################################

for sample in `ls ${pathHisat}/ `
do
    if [ -e ${pathResults}/${sample}/CoverageVariation_Ensembl${release}_20windows_minlength400.txt ]; then 
	echo "already done"
    else
	if [ -e ${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz ]&&[ -e ${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz ]; then
	    echo "doing "${sample}
	    
	    if [ -e ${pathResults}/${sample} ]; then
		echo "path exists"
	    else
		mkdir ${pathResults}/${sample}
	    fi
	    
	    echo "#!/bin/bash" > bsub_script_covvar
	  
	    echo "perl ${pathScripts}/coverage.variation.pl --pathExonBlocks=${pathEnsembl}/ExonBlocks_FilteredTranscripts_Ensembl${release}.txt --pathCoverageForward=${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz  --pathCoverageReverse=${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz --nbWindows=20 --minLength=400 --pathOutput=${pathResults}/${sample}/CoverageVariation_Ensembl${release}_20windows_minlength400.txt"  >> bsub_script_covvar
	    
	    qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_covvar_${sp}_${sample}.txt -e ${pathScripts}/std_error_${sp}_${sample}.txt ${pathScripts}/bsub_script_covvar

	fi
    fi
done

######################################################################################################
