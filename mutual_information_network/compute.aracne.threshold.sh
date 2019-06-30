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

java -Xmx4G -jar ${pathAracne}/dist/aracne.jar -e ${pathResults}/${sp}/TPM.txt --tfs ${pathResults}/${regulators} -o ${pathResults}/${sp}/ --pvalue 1E-8 --seed 1 --calculateThreshold

####################################################################################
