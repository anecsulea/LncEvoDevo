#!/bin/bash

################################################################

export path=LncEvoDevo
export pathAln=${path}/data/genome_alignments
export pathScripts=${path}/scripts/coding_potential/CSF

################################################################

export chrlist=""

for i in {1..20} X Y
do
    export chrlist=chr${i},${chrlist}
done

################################################################

perl ${pathScripts}/divide.maf.pl --pathAlignment=${pathAln}/Rat_20way/rn6.20way.maf.gz --refSpecies=rn6 --chromosomes=${chrlist} --dirOutput=${pathAln}/Rat_20way/

################################################################
