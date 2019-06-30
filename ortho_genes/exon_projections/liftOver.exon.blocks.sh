#!/bin/bash

export ref=$1
export tg=$2
export annot=$3

######################################################################


export path=LncEvoDevo
    
######################################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathStringTie=${path}/results/stringtie_assembly
export pathGenomeAlignments=${path}/data/genome_alignments
export pathResults=${path}/results/exon_projections
export pathScripts=${path}/scripts/ortho_genes/exon_projections

export release=94

######################################################################

export pathLiftOver=""

if [ ${ref} = "Mouse" ]&&[ ${tg} = "Rat" ]; then
    export pathLiftOver=${pathGenomeAlignments}/mm10ToRn6.over.chain.gz  
fi

if [ ${ref} = "Rat" ]&&[ ${tg} = "Mouse" ]; then
    export pathLiftOver=${pathGenomeAlignments}/rn6ToMm10.over.chain.gz  
fi

if [ ${ref} = "Mouse" ]&&[ ${tg} = "Chicken" ]; then
    export pathLiftOver=${pathGenomeAlignments}/mm10ToGalGal5.over.chain.gz  
fi

if [ ${ref} = "Chicken" ]&&[ ${tg} = "Mouse" ]; then
    export pathLiftOver=${pathGenomeAlignments}/galGal5ToMm10.over.chain.gz  
fi

if [ ${ref} = "Rat" ]&&[ ${tg} = "Chicken" ]; then
    export pathLiftOver=${pathGenomeAlignments}/rn6ToGalGal5.over.chain.gz  
fi

if [ ${ref} = "Chicken" ]&&[ ${tg} = "Rat" ]; then
    export pathLiftOver=${pathGenomeAlignments}/galGal5ToRn6.over.chain.gz  
fi


if [ ${ref} = "Mouse" ]&&[ ${tg} = "Human" ]; then
    export pathLiftOver=${pathGenomeAlignments}/mm10ToHg38.over.chain.gz  
fi

if [ ${ref} = "Rat" ]&&[ ${tg} = "Human" ]; then
    export pathLiftOver=${pathGenomeAlignments}/rn6ToHg38.over.chain.gz  
fi


if [ ${pathLiftOver} = "" ]; then  
    echo "don't know what chain to use for "${ref} ${tg}
    exit
fi

######################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}/${ref}/bed_parts
    export prefixExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}/${ref}/combined/bed_parts
    export prefixExons=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
fi

######################################################################

for file in `ls ${pathExons} | grep ${prefixExons} | grep part`
do
    export prefix=`basename ${file} .bed`
    echo ${file} ${prefix}
  
    if [ -e  ${pathResults}/From${ref}_To${tg}_${prefix}.bed ]; then
	echo "already done"
    else
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_liftOver
	
	
	echo "liftOver -minMatch=0.1 -multiple ${pathExons}/${prefix}.bed ${pathLiftOver} ${pathResults}/From${ref}_To${tg}_${prefix}.bed ${pathResults}/From${ref}_To${tg}_${prefix}.unmapped " >> ${pathScripts}/bsub_script_liftOver


	qsub -q long -l s_rss=4G,sps=1 -o ${pathScripts}/std_output_liftOver_${ref}_${tg}_${annot}.txt -e ${pathScripts}/std_error_liftOver_${ref}_${tg}_${annot}.txt ${pathScripts}/bsub_script_liftOver
	
    fi
done 

######################################################################
