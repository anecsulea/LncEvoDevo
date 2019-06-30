#!/bin/bash

export species=$1
export annot=$2
export phast=$3

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathResults=${path}/results/sequence_evolution/phastcons/${species}/${phast}
export pathPhastCons=${path}/data/phastcons/${species}
export pathScripts=${path}/scripts/sequence_evolution/phastcons

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}
    export suffixExons=IntronBlocks_AllGenes_FilteredTranscripts_Ensembl${release}
    export suffixActualExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}
    export suffixExons=IntronBlocks_AllGenes_FilteredTranscripts_StringTie_Ensembl${release}
    export suffixActualExons=ExonBlocks_FilteredTranscripts_StringTie_Ensembl${release}
fi

#####################################################################

if [ -e ${pathResults} ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}
fi

#####################################################################

if [ ${species} = "Mouse" ]; then
    if [ ${phast} = "placental" ]; then
	export suffix=.phastCons60way.placental.wigFix.gz 
    fi

    if [ ${phast} = "60vertebrates" ]; then
	export suffix=.phastCons60way.wigFix.gz 
    fi
fi

#####################################################################

for chr in {1..19} X Y 
do
    export pathPhast=${pathPhastCons}/chr${chr}${suffix}

    if [ -e ${pathPhast} ]; then
    
	echo "#!/bin/bash" > bsub_script


	echo "perl ${pathScripts}/compute.phastcons.exons.pl --pathExonBlocks=${pathExons}/${suffixExons}.txt --pathMaskExonBlocks=${pathExons}/${suffixActualExons}.txt --pathPhastCons=${pathPhast} --chr=${chr} --pathOutputExons=${pathResults}/PhastCons_IntronAverage_MaskedExons_chr${chr}_${suffixExons}.txt --pathOutputGenes=${pathResults}/PhastCons_GeneAverage_MaskExons_chr${chr}_${suffixExons}.txt " >> bsub_script

	
	if [ ${chr} = Y ]||[ ${chr} = 19 ]||[ ${chr} = 20 ]; then
	    qsub -P P_biometr -q huge -l s_rss=10G,sps=1  -o ${pathScripts}/std_output_phastcons_introns_noov_${species}.txt -e ${pathScripts}/std_error_phastcons_introns_noov_${species}.txt ${pathScripts}/bsub_script	
	else
	    qsub -P P_biometr -q mc_highmem_long -l s_rss=30G,sps=1 -pe multicores 1 -o ${pathScripts}/std_output_phastcons_introns_noov_${species}.txt -e ${pathScripts}/std_error_phastcons_introns_noov_${species}.txt ${pathScripts}/bsub_script	
	fi
	
    fi
done

#####################################################################

   
