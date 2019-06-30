#!/bin/bash

export sp=$1
export alntype=$2
export windowSize=$3

##################################################################

export path=LncEvoDevo
export pathMatrices=${path}/results/CSF/matrices/${alntype}/${sp}
export pathResults=${path}/results/CSF/genome_scores/${alntype}/${sp}
export pathScripts=${path}/scripts/coding_potential/CSF
export pathGenomeAlignments=${path}/data/genome_alignments/${alntype}
export pathGeneticCode=${pathScripts}/standard_genetic_code.txt 

##################################################################

export penalty=0.001
export pseudofreq=0.00001

export maxsize=300000000

if [ ${sp} = "Opossum" ]; then
    export maxsize=800000000
fi 

echo ${sp}" max size "${maxsize}

##################################################################

export paths=""

for chr in {1..32} X Y Z W
do
    if [ -e ${pathGenomeAlignments}/chr${chr}.maf.gz ]; then
	
	export start="0"
	
	while [ ${start} -lt ${maxsize} ]
	do
	    export end=$[${start}+20000000]

	    if [ -e ${pathResults}/CSF_CoveredRegions_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt ]; then
		export paths=${pathResults}/CSF_CoveredRegions_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}_chr${chr}_start${start}_end${end}.txt,${paths}
	    fi
	    
	    export start=$[${start}+20000000]
	done
    fi
done

##################################################################

perl ${pathScripts}/combine.covered.regions.pl --pathsRegions=${paths} --pathOutput=${pathResults}/CSF_CoveredRegions_windowSize${windowSize}bp_penalty${penalty}_pseudofreq${pseudofreq}.txt 

##################################################################
