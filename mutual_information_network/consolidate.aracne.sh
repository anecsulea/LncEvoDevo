#!/bin/bash

export sp=$1
export type=$2

####################################################################################

export path=LncEvoDevo
export pathTools=Tools

export pathResults=${path}/results/mutual_information_network/${type}
export pathAracne=${pathTools}/ARACNe-AP/
export pathScripts=${path}/scripts/mutual_information_network

####################################################################################

if [ ${type} = "all_genes" ]; then
    export regulators=AllGenes.txt
fi

if [ ${type} = "lncRNAs_only" ]; then
    export regulators=LncRNAs.txt
fi

####################################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_consolidate

echo "java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -o ${pathResults}/${sp}/ --consolidate --nobonferroni" >>  ${pathScripts}/bsub_script_consolidate

# qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_aracne_${sp}.txt -e ${pathScripts}/std_error_aracne_${sp}.txt ${pathScripts}/bsub_script_consolidate

qsub -q mc_highmem_huge -l s_rss=30G,sps=1 -o ${pathScripts}/std_output_aracne_${sp}.txt -e ${pathScripts}/std_error_aracne_${sp}.txt ${pathScripts}/bsub_script_consolidate

####################################################################################
