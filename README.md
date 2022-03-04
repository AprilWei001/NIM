# NIM

## Table of contents
* [Whole-genome simulation](#Whole-genome-simulation)
* [SNP annotation](#SNP-annotation)
* [Estimating heritability with RHE-mc](#Estimating-heritability-with-RHE-mc)
* [Fine mapping](#Fine-mapping)


## Whole-genome simulation
Simulation cript: `simAnyArchitecture.sh`
*This simulation script depends on software: [plink](https://www.cog-genomics.org/plink2/), [gcta64](https://cnsgenomics.com/software/gcta/#Overview), script: `shuffle.py`*

*input files: <genotype input file name (e.g. qced)>  <.frq> which has the in-sample MAF of SNPs used, and <.ld>*

In-sample ldscore is comptued with [gcta64](https://cnsgenomics.com/software/gcta/#Overview) using all QC-ed genotypes

      gcta64 --bfile <genotype input file name (e.g. qced)> --ld-score --ld-wind 10000 --out <output filename> 

With all dependencies in place, run simulation with:

    ./simAllArchitecture.sh <type of architecture> <seed>
for example 

    ./simAllArchitecture.sh POLY 1 
will give the simulation for polygenic architecture with seed 1. 

The types of architecture can be:
* POLY -- polygenic architecture 
* LOW -- low LD SNPs (defined as ldsc <=10) contributes to 9000 causal variants, and high LD SNPs contributes 1000 causal variants, 
* HIGH -- high LD (defined as ldsc > 10) contributes to 9000 causal variants, and low LD SNPs contributes 1000 causal variants,
* RARE -- rare SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and common SNPs contributes 1000 causal variants, 
* COMMON -- common SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and rare SNPs contributes 1000 causal variants,
* 
Within simAllArchitecture.sh
hsq="0.5" or "0.2" for high and low heritability simulation
ncausal=10000 or 100000 for moderate and high polygenecity simulation

## SNP annotation
### Ancestry 
We use two annotation of ancestry: neanderthal ancestry and modern human ancestry. 
We first get the list of confident NIMs from [Sankararaman et al 2014](https://www.nature.com/articles/nature12961) in 1KG EUR populations. Then we expand it to include all the SNPs in strong LD (r^2 >= 0.99) with any confident NIMs as expanded NIMs with `plink --bfile <qced>  --show-tags <confident NIMs>  --tag-r2 0.99  --tag-kb 200 --out <expanded NIMs>`
All the SNPs that are not expanded NIMs are annotated as modern human ancestry.

### MAF
We use 5 MAF bins, which split all QC-ed SNPs from low to high MAF into equal sized bins
### LD
We use 5 LD bins, which split all QC-ed SNPs from low to high LD-score into equal sized bins.

### Consruct non-overlapping annotation

In the output file name annotation for ancestry is indicated as '.anc', maf with '.maf', and ld-score with '.ld'. 

The annotation is non-overlapping, meaning that SNPs are partitioned into combinations of annotation, hence each SNP belongs to one and only one annotation. 

For example, in a file named with 'nol.anc.maf.annot', there are 10 annotations total (2 ancestry * 5 maf), with the Neanderthal ancestry into 5 maf bins, and modern human ancestry into 5 maf bins. 

*An exception is '.anc.maf.ld' annotation which should have 2*5*5 = 50 non-overlapping bins, but because there are few high MAF Neanderthal SNPs and low LD-score Neanderthal SNPs, some bins are empty or near empty, hence we remove the bins with smaller than 30 SNPs in them. This particular annotation is only constructed for Tag SNPs because the number of confident NIMs is too small to have so many annotations.*

## Estimating heritability with RHE-mc
### Information about RHE-mc
In this study we use an extension of RHE-mc which takes a coefficient file to allow us to define new summary statistics from linear combinations (`-coeff flag`) of the variance components. Both the point estimate and the standard error of the summary statistics are estimated from jackknife.
Information about the original RHE-mc can be found at [RHE-mc GitHub](https://github.com/sriramlab/RHE-mc)
This new version with extension is now available at [RHE-mc add link](https:)
### Simulated data
Note that, the gcta64 simulated phenotype does not have a header, e.g.:

        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 

This file (without header) can be directly used with plink for GWAS, but RHE-mc takes phenotype file which starts with a header, hence require adding a header line with e.g.:

        sed  -i '1 i\FID IID pheno' *.phen
Then we can get the correct format for RHE-mc, e.g.:

        FID IID pheno
        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 

Run RHE-mc with the supply of genotype, phenotype, and annotation files for whole-genome simulated data with

        RHEmc_mem -g <genotype file> -p <phenotype file> -coeff <weight coefficient> -annot <annotation file> -k 10 -jn 100  -o <output file>

### UKBB data       
Run RHE-mc with the supply of genotype, phenotype, coavariate and anonotation files for UKBB data.

        RHEmc_mem -g <genotype file>  -p <phenotype file> -c <covar file>  -coeff <weight coefficient>  -annot <annotation file> -k 10 -jn 100  -o <output file>
   
### H2 Partitioning
From the output of RHE-mc, we extract the heritabiliity and standard error of heritability for each summary statatistic we defined. The output starts with the coefficients we defined, followed by the statistics we defined, and more, for example:

        Coefficients statistics  :
          1           1           1           1           1           1           1           1           1           1
          1           1           1           1           1           0           0           0           0           0
          0           0           0           0           0  0.00293944  0.00971893   0.0229379   0.0101663 0.000217432
         -1          -1          -1          -1          -1  0.00293944  0.00971893   0.0229379   0.0101663 0.000217432

        0-th statistic:  point estimate: 0.189033 SE: 0.00626612
        1-th statistic:  point estimate: 0.00109391 SE: 0.000826601
        2-th statistic:  point estimate: 0.0020683 SE: 0.000117277
        3-th statistic:  point estimate: 0.000974397 SE: 0.000885459
        OUTPUT: 
        Variances: 
        Sigma^2_0: 7.73175e-05 ,  se: 0.000205804
        Sigma^2_1: -0.000235195 ,  se: 0.00029295
        ... ...
        ... ...

## Fine mapping

### Information about SuSiE 
Information about Susie can be found at [SuSiE GitHub](https://stephenslab.github.io/susie-paper/index.html). It is an R package for fine mapping.
        
        
### Benchmark SuSiE with simulated data

1. GWAS with plink, 

        plink --silent --bfile <qced>--pheno $OUT.phen --linear standard-beta --out <output file>
        
  then extract significant SNPs with P < 10^(-10) with `getSignificantSNPs.m`
  
2. LD pruning of significant SNPs 

        plink --bfile <qced>  --extract <list of SNPs>  --indep-pairwise 100kb 1 0.99 --out <output filename>
        
3. For each pruned-in GWAS signficant SNP
    3.1 Output the 200 kb region surrounding this SNP, for example:

        plink --bfile <qced> --snps <range of the snps (e.g. 1:4669624-1:4869948)> --recode A --out <output genotype text file>
        
    3.2 Run SuSiE on the 200kb region surrounding the SNP with script `runSusie.R` 

    *`runSusie.R` takes the genotype text file input from the tested region generated in 3.1*

        R --slave --args <genotype text file > <phenotype file> < runSusie.R >
        
    *SuSiE requres the genotype file to be loaded in R, hence asks for quite a bit of memory. For 200kb region in the UKBB data, it could be as small as 20 GB or as large as 50 GB.*
    3.3 Remove the genotype data and plink log file
    
*step 3 is automated by first outputing the range of all pruned-in SNP (with `getSuisieRange.m`) then looping through 3.1-3.3.* with bash script
4. Analyze SuSiE outputs and determine the boundary and credible NIM set

### SuSiE to UKBB phenotypes
The only difference between handling UKBB phenotypes and simulated phenotypes is that there is covariates in UKBB. So we first compute the phenotypic residuals by linear regression to regress out the covariates (e.g., PC, NIM PC, sex, age), and then normalize these residuals, before using them as input to SuSiE at step 3.2.
