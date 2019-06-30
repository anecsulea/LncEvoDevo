#!/bin/bash

export sp=$1
export aln=$2
export annot=$3
export winsize=$4

##############################################################

export path=LncEvoDevo

##############################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathCSF=${path}/results/CSF/genome_scores/${aln}/${sp}
export pathRepeats=${path}/results/overlap_repeats/${sp}
export pathResults=${path}/results/coding_potential/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF

export release=94

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}/FilteredTranscripts_Ensembl${release}.gtf
    export prefixOutput=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export pathGTF=${pathStringTie}/FilteredTranscripts_StringTie_Ensembl${release}.gtf
    export prefixOutput=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_class

echo "perl ${pathScripts}/classify.genes.allexons.pl --pathAnnotGTF=${pathGTF} --pathPositiveScores=${pathCSF}/CSF_PositiveScores_windowSize${winsize}bp_penalty0.001_pseudofreq0.00001.txt --pathCoveredRegions=${pathCSF}/CSF_CoveredRegions_windowSize${winsize}bp_penalty0.001_pseudofreq0.00001.txt --minFractionOverlap=0.0 --minLengthOverlap=150 --pathOutputExons=${pathResults}/OverlapExons_AllExons_CSFScores_${aln}_windowSize${winsize}_${prefixOutput}.txt --pathOutputTranscripts=${pathResults}/TranscriptClassification_AllExons_CSFScores_${aln}_windowSize${winsize}_${prefixOutput}.txt  --pathOutputGenes=${pathResults}/GeneClassification_AllExons_CSFScores_${aln}_windowSize${winsize}_${prefixOutput}.txt" >>  ${pathScripts}/bsub_script_class

qsub -q huge -l s_rss=6G,sps=1 -o ${pathScripts}/std_output_classCSF_${sp}_${aln}_${annot}_${winsize}.txt -e ${pathScripts}/std_error_classCSF_${sp}_${aln}_${annot}_${winsize}.txt ${pathScripts}/bsub_script_class


##############################################################
