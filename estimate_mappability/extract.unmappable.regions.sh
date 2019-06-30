#/bin/bash

export species=$1

##############################################################

export path=LncEvoDevo
export pathResults=${path}/results/mappability/${species}
export pathScripts=${path}/scripts/estimate_mappability

##############################################################

export pathsMap=""

for chr in {1..32} X Y Z W
do
    
    if [ -e ${pathResults}/fake_reads_chr${chr}_mappedregions.txt ]; then
	export pathsMap=${pathResults}/fake_reads_chr${chr}_mappedregions.txt,${pathsMap}
    else
	if [ -e ${pathResults}/fake_reads_chr${chr}.sam.gz ]; then
	    echo "cannot find "${pathResults}/fake_reads_chr${chr}_mappedregions.txt
	    exit
	fi
    fi
done

##############################################################

perl ${pathScripts}/extract.unmappable.regions.pl --pathsMappedRegions=${pathsMap} --pathChromosomeSizes=${pathResults}/chromosome_sizes_samtools.txt --pathOutputTxt=${pathResults}/unmappable_regions.txt --pathOutputBedGraph=${pathResults}/unmappable_regions.bedGraph

##############################################################
