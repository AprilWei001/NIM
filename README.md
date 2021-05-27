# NIM

## Table of contents
* [UKBB SNP QC](#UKBB-SNP-QC)
* [Whole-genome simulation](#Whole-genome-simulation)
* [Estimating heritability with RHE-mc](#Estimating-heritability-with-RHE-mc)
* [Fine mapping](#Fine-mapping)


## UKBB SNP QC 
SNP QC for UKBB whole-genome imputed data with plink flags `--exclude <SNPs in the mhc region txt file> --maf 0.001 --geno 0.01 --hwe 1e-7 --keep <ukbb wb unrelated individuals fam file> --make-bed --out <output file name (e.g., qced)>`

%{note to self: find this information on hoffman /u/project/sgss/UKBB/data/imp/Qced/All/README, although the QC-ed data used are in /u/project/sgss/UKBB/data/imp/Qced/maf_0.001 which persumably only changes the flag --maf 0.01 into --maf 0.001%}

## Whole-genome simulation
Script: simAnyArchitecture.sh
* This simulation depends on software: plink, gcta64, cripts: shuffle.py
* input files: <genotype input file name (e.g. qced)>  <.frq> which has the in-sample MAF of SNPs used, and <.ld> in-sample ldscore: comptued with `gcta64 --bfile <genotype input file name (e.g. qced)> --ld-score --ld-wind 10000 --out <output filename>` 

Run with: `./simAllArchitecture.sh <type of architecture> <seed>`
  for example `./simAllArchitecture.sh POLY 1` will give the simulation for polygenic architecture with seed 1. 

The types of architecture can be:
* POLY -- polygenic architecture 
* LOW -- low LD SNPs (defined as ldsc <=10) contributes to 9000 causal variants, and high LD SNPs contributes 1000 causal variants, 
* HIGH -- high LD (defined as ldsc > 10) contributes to 9000 causal variants, and low LD SNPs contributes 1000 causal variants,
* RARE -- rare SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and common SNPs contributes 1000 causal variants, 
* COMMON -- common SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and rare SNPs contributes 1000 causal variants,

%{note to self: script and small sized dependencies moved to a folder on tabla: /home/aprilwei/projects/nimHeretability/github -- original simulations are outputs from other scripts but merged here for the purpose of sharing, tested runnable on tabla.%}

## Estimating heritability with RHE-mc

### SNP annotation

### Partitioning

### META-analysis

## Fine mapping
