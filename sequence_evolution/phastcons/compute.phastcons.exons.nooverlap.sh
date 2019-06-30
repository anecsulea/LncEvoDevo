#!/bin/bash

export species=$1
export phast=$2

#####################################################################

export path=LncEvoDevo

export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathGeneOverlaps=${path}/results/gene_overlaps/${species}
export pathPhastCons=${path}/data/phastcons/${species}
export pathResults=${path}/results/sequence_evolution/phastcons/${species}/${phast}
export pathScripts=${path}/scripts/sequence_evolution/phastcons

export release=94

#####################################################################

export pathExons=${pathGeneOverlaps}
export suffixExons=ExonBlocks_ExcludingOverlapOtherGenes_FilteredTranscripts_StringTie_Ensembl${release}

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

	
	echo "perl ${pathScripts}/compute.phastcons.exons.pl --pathExonBlocks=${pathExons}/${suffixExons}.txt --pathPhastCons=${pathPhast} --chr=${chr} --pathOutputExons=${pathResults}/PhastCons_ExonAverage_chr${chr}_${suffixExons}.txt --pathOutputGenes=${pathResults}/PhastCons_GeneAverage_chr${chr}_${suffixExons}.txt " >> bsub_script
	
	if [ ${chr} = Y ]||[ ${chr} = 19 ]||[ ${chr} = 20 ]; then
	    qsub -P P_biometr -q huge -l s_rss=10G,sps=1  -o ${pathScripts}/std_output_phastcons_exons_noov_${species}.txt -e ${pathScripts}/std_error_phastcons_exons_noov_${species}.txt ${pathScripts}/bsub_script	
	else
	    qsub -P P_biometr -q mc_highmem_long -l s_rss=30G,sps=1 -pe multicores 1 -o ${pathScripts}/std_output_phastcons_exons_noov_${species}.txt -e ${pathScripts}/std_error_phastcons_exons_noov_${species}.txt ${pathScripts}/bsub_script	
	fi
	
    fi
done


#####################################################################

