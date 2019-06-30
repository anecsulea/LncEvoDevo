#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3

####################################################################################

export path=LncEvoDevo

####################################################################################

export pathLiftOver=${path}/data/genome_alignments
export pathProjections=${path}/results/exon_projections
export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathUCSC=${path}/data/UCSC_sequences
export pathResults=${path}/results/ortho_genes/whole_genome_alignments
export pathScripts=${path}/scripts/ortho_genes/whole_genome_alignments

export release=94

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
    export suffix=""
fi

if [ ${annot} = "StringTie" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
    export pathExons=${pathStringTie}
    export suffix="/combined"
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

####################################################################################

if [ -e ${pathResults}/tba_alignments_projection_clusters/${sp1}_${sp2}_${annot} ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}/tba_alignments_projection_clusters/${sp1}_${sp2}_${annot}
fi

####################################################################################

for i in {0..9}
do    
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_tba
    
    echo "perl ${pathScripts}/run.tba.alignments.projection.clusters.pl --species1=${sp1} --species2=${sp2} --pathSequences1=${pathExons}/${sp1}${suffix}/${prefix}_cDNASequences.fa --pathSequences2=${pathExons}/${sp2}${suffix}/${prefix}_cDNASequences.fa --pathClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefix}.txt --run=${i} --modulo=10 --dirOutput=${pathResults}/tba_alignments_projection_clusters/${sp1}_${sp2}_${annot}" >> ${pathScripts}/bsub_script_tba
    
    qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_tba_${sp1}_${sp2}_${annot}_${i}.txt -e ${pathScripts}/std_error_tba_${sp1}_${sp2}_${annot}_${i}.txt ${pathScripts}/bsub_script_tba
    
done
####################################################################################
