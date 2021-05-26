# NIM

## Table of contents
* [UKBB SNP QC](#UKBB-SNP-QC)
* [Whole genome simulation](#Whole-genome-simulation)


## UKBB SNP QC 
QC with plink flags `--exclude <SNPs in mhc region txt file> --maf 0.001 --geno 0.01 --hwe 1e-7 --keep <ukbb wb unrelated individuals fam file> --make-bed --out <output file name (e.g., qced)>`

## Whole genome simulation
Simulate any genetic architecture on tabla: /home/aprilwei/projects/nimHeretability/wgSim/simAllArchitecture.sh 

This simulation depends on software: plink, gcta64, cripts: shuffle.py
input files: <genotype input file name (e.g. qced)>  <.frq> which has the in-sample MAF of SNPs used, and <.ld> in-sample ldscore: comptued with `gcta64 --bfile <genotype input file name (e.g. qced)> --ld-score --ld-wind 10000 --out <output filename>` 

Run with `./simAllArchitecture.sh <type of architecture> <seed>` for example `./simAllArchitecture.sh POLY 1` will give the simulation for polygenic architecture with seed 1. The type of architecture can be "POLY", "LOW", "HIGH", "RARE", "COMMON".

