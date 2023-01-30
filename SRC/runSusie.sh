#!/bin/bash



#!/bin/bash


for pcutoff in 10; do
for pheno in `cat filename96pheno.txt`; do        
	filename="susieRange/$pheno.pcutoff.$pcutoff.signim.nimrange";
	if test -f "$filename"; then       	
		num=`wc -l $filename |awk '{print $1}'`;
		for index in `seq 1 $num`; do
			out=susie/$pheno.pcutoff.$pcutoff.$index.susie
			raw=$pheno.pcutoff.$pcutoff.$index.raw
			if test ! -f "$out"; then
				echo "./printRGenotypeUKBB.sh $pcutoff $pheno $index $filename";
			elif test -f "$raw"; then
				echo "./printRGenotypeUKBB.sh $pcutoff $pheno $index $filename";	
			fi 
      		done
	fi
done |msub -R h_data=50G,h_rt=10:00:00,highp -j susieUKBB.$pcutoff -g 20;
done
