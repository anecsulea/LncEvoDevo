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
export pathGeneOverlaps=${path}/results/gene_overlaps/${sp}
export pathResults=${path}/results/coding_potential/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF

export release=94

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export prefixAnnot=FilteredTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "StringTie" ]; then
    export prefixAnnot=FilteredTranscripts_StringTie_Ensembl${release}
fi

##############################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_class

echo "perl ${pathScripts}/classify.genes.nonoverlap.pl --pathNonOverlappingAnnot=${pathGeneOverlaps}/ExonCoordinates_ExcludingOverlapOtherGenesRepeatsRetrogenes_${prefixAnnot}.txt --pathPositiveScores=${pathCSF}/CSF_PositiveScores_windowSize${winsize}bp_penalty0.001_pseudofreq0.00001.txt --pathCoveredRegions=${pathCSF}/CSF_CoveredRegions_windowSize${winsize}bp_penalty0.001_pseudofreq0.00001.txt --minFractionOverlap=0.0 --minLengthOverlap=150 --pathOutputExons=${pathResults}/OverlapExons_NonOverlappingRegions_CSFScores_${aln}_windowSize${winsize}_${prefixAnnot}.txt --pathOutputTranscripts=${pathResults}/TranscriptClassification_NonOverlappingRegions_CSFScores_${aln}_windowSize${winsize}_${prefixAnnot}.txt  --pathOutputGenes=${pathResults}/GeneClassification_NonOverlappingRegions_CSFScores_${aln}_windowSize${winsize}_${prefixAnnot}.txt" >>  ${pathScripts}/bsub_script_class

qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_classCSF_${sp}_${aln}_${annot}_${winsize}.txt -e ${pathScripts}/std_error_classCSF_${sp}_${aln}_${annot}_${winsize}.txt ${pathScripts}/bsub_script_class

##############################################################
