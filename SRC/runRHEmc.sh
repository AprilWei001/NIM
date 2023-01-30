#!/bin/bash

pathGeno="../all";
annot="expanded.anc";
pathAnnot="$annot";


for item in `cat filename96pheno.txt`; do 
echo "RHEmc_mem -g $pathGeno  -p ../ukbbPheno/$item -c ukbb.nim.covar -annot $pathAnnot.annot -coeff $pathAnnot.weight -k 10 -jn 100  -o rhemc/$item.$annot.rhemc > logs/$item.$annot.log";
done |msub -R h_data=24G,h_rt=4:00:00:00,highp -j $annot -g 1;

