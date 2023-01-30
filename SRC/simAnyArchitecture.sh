#!/bin/bash

TYPE=$1
SEED=$2
hsq=$3
ncausal=$4

path="/home/data/ukbiobank/data/geno/imp/hard_calls_291273/" #local path to the genomic location
#hsq="0.2"   #for wg, use 0.5
#ncausal=100000
n0=`echo $ncausal|awk '{printf "%d\n" ,$1/10}'`
n1=`echo $ncausal|awk '{printf "%d\n" ,9*$1/10}'`

OUT="ncausal.$ncausal.hsq.$hsq.$TYPE.$SEED"

if [ "$TYPE" == "BASELINE" ]; then
cat $path/qced.bim | awk '{ print $2 }' | python ../shuffle.py $SEED | head -n $ncausal > $OUT.causal
elif [ "$TYPE" == "RARE" ]; then
cat qced.frq | awk '$5 <= 0.05 { print $2 }' | python ../shuffle.py $SEED | head -n $n1 > $OUT.causal
cat qced.frq | awk '$5 > 0.05 { print $2 }' | python ../shuffle.py $SEED | head -n $n0 >> $OUT.causal
elif [ "$TYPE" == "COMMON" ]; then
cat qced.frq | awk '$5 <= 0.05 { print $2 }' | python ../shuffle.py $SEED | head -n  $n0 > $OUT.causal
cat qced.frq | awk '$5 > 0.05 { print $2 }' | python ../shuffle.py $SEED | head -n $n1 >> $OUT.causal
elif [ "$TYPE" == "LOW" ]; then
cat all.score.ld | awk '$8 <= 10 { print $1 }' | python ../shuffle.py $SEED | head -n $n1 > $OUT.causal
cat all.score.ld | awk '$8 > 10 { print $1 }' | python ../shuffle.py $SEED | head -n $n0 >> $OUT.causal
elif [ "$TYPE" == "HIGH" ]; then
cat all.score.ld | awk '$8 <= 10 { print $1 }' | python ../shuffle.py $SEED | head -n  $n0 > $OUT.causal
cat all.score.ld | awk '$8 > 10 { print $1 }' | python ../shuffle.py $SEED | head -n $n1 >> $OUT.causal
fi

/home/sriram/bin/plink --bfile $path/qced --extract $OUT.causal --make-bed --out $OUT

/home/sriram/bin/gcta/gcta64 --bfile $OUT --simu-qt --simu-hsq $hsq --out $OUT --extract $OUT.causal --simu-causal-loci $OUT.causal #edited

# COMPUTE ASSOCIATION STATISTICS
/home/sriram/bin/plink --silent --bfile $path/qced \
--pheno $OUT.phen --linear standard-beta --out $OUT
