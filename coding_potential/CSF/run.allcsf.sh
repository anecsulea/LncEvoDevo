#!/bin/bash

##############################################################

for annot in Ensembl StringTie
do
    ./classify.genes.nonoverlap.sh Mouse Mouse_60way ${annot} 75 local 
    ./classify.genes.allexons.sh Mouse Mouse_60way ${annot} 75 local
done 

for annot in Ensembl StringTie
do
    ./classify.genes.nonoverlap.sh Rat Human_100way ${annot} 75 local
    ./classify.genes.allexons.sh Rat Human_100way ${annot} 75 local
done 


for annot in Ensembl StringTie
do
    ./classify.genes.nonoverlap.sh Rat Rat_20way ${annot} 75 local
    ./classify.genes.allexons.sh Rat Rat_20way ${annot} 75 local
done 


for annot in Ensembl StringTie
do 
    ./classify.genes.nonoverlap.sh Chicken Human_100way ${annot} 75 local
    ./classify.genes.allexons.sh Chicken Human_100way ${annot} 75 local
done 

##############################################################
