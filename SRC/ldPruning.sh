#!/bin/bash
for pheno in `cat sigNIM.txt`; do
if test ! -f "prunedNIM/$pheno.prune.in";then
echo "~/bin/./plink --bfile ../nim --extract sigNIM/$pheno --indep-pairwise 100kb 1 0.99 --out prunedNIM/$pheno";
fi
done |msub  -R h_data=12G,h_rt=2:00:00 -j ldpruning -g 4
