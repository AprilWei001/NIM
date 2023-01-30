#!/bin/bash
path="/u/project/sriram/crrobles/ukbb/results/imp.ukb0.5/april/nims_pccorr";

pcutoff=10
for pheno in `cat filename96pheno.txt`; do
cat $path/$pheno.*glm* |awk '$12<1E-10 {print $3}' > sigNIM/$pheno.pcutoff.$pcutoff.signim
done


