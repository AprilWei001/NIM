#this is the snpEFF script to make the Data S8 file.  A VCF for all the credible nims must be created first.
java -Xmx8g -jar snpEff.jar -o bed -v -csvStats ex1.csv GRCh37.75 ./all.pheno.all.credible.covar.vcf > Data\ S8.txt
