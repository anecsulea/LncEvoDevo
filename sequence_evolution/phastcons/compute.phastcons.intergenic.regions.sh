#!/bin/bash

export species=$1
export annot=$2
export phast=$3

#####################################################################

export path=LncEvoDevo

export pathResults=${path}/results/sequence_evolution/phastcons/${species}/${phast}
export pathStringTie=${path}/results/stringtie_assembly/${species}/combined
export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathPhastCons=${path}/data/phastcons/${species}
export pathScripts=${path}/scripts/sequence_evolution/phastcons

export release=94

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathExons=${pathEnsembl}
    export suffixExons=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "StringTie" ]; then
    export pathExons=${pathStringTie}
    export suffixExons=FilteredTranscripts_StringTie_Ensembl${release}
 fi

export minsize=1kb
export mindist=5kb

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

for chr in  {1..19} X Y
do
    export pathPhast=${pathPhastCons}/chr${chr}${suffix}

    if [ -e ${pathPhast} ]; then

	if [ -e ${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_chr${chr}_${suffixExons}.txt ]; then
	    echo "already done"
	else
    	    echo "#!/bin/bash" > bsub_script
	 
	    
	    echo "perl ${pathScripts}/compute.phastcons.intergenic.regions.pl --pathIntergenicRegionCoords=${pathExons}/IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_${suffixExons}.txt --pathPhastCons=${pathPhast} --chr=chr${chr} --pathOutput=${pathResults}/PhastCons_IntergenicRegions_MinDistance${mindist}_MinSize${minsize}_chr${chr}_${suffixExons}.txt" >> bsub_script

	    if [ ${chr} = Y ]||[ ${chr} = 19 ]||[ ${chr} = 20 ]; then
		qsub -P P_biometr -q huge -l s_rss=10G,sps=1  -o ${pathScripts}/std_output_phastcons_promoters_${species}.txt -e ${pathScripts}/std_error_phastcons_promoters_${species}.txt ${pathScripts}/bsub_script	
	    else
		qsub -P P_biometr -q mc_highmem_long -l s_rss=30G,sps=1 -pe multicores 1 -o ${pathScripts}/std_output_phastcons_promoters_${species}.txt -e ${pathScripts}/std_error_phastcons_promoters_${species}.txt ${pathScripts}/bsub_script	
	    fi
	    	   
	fi
    fi
done

#####################################################################

   
