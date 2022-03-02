# NIM

## Table of contents
* [UKBB SNP QC](#UKBB-SNP-QC)
* [Whole-genome simulation](#Whole-genome-simulation)
* [SNP annotation](#SNP-annotation)
* [Estimating heritability with RHE-mc](#Estimating-heritability-with-RHE-mc)
* [Fine mapping](#Fine-mapping)


## UKBB SNP QC 
SNP QC for UKBB whole-genome imputed data is the same as described in [RHE-mc paper](https://www.nature.com/articles/s41467-020-17576-9)

        plink  --bfile <ukbb imputed genotype data> --exclude <SNPs in the mhc region txt file> --maf 0.001 --geno 0.01 --hwe 1e-7 --keep <ukbb wb unrelated individuals fam file> --make-bed --out <output file name (e.g., qced)>


## Whole-genome simulation
Simulation cript: `simAnyArchitecture.sh`
*This simulation script depends on software: [plink](https://www.cog-genomics.org/plink2/), [gcta64](https://cnsgenomics.com/software/gcta/#Overview), script: `shuffle.py`*

*input files: <genotype input file name (e.g. qced)>  <.frq> which has the in-sample MAF of SNPs used, and <.ld>*

In-sample ldscore is comptued with 

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



## SNP annotation
### Ancestry 
We use two annotation of ancestry: neanderthal ancestry and modern human ancestry. We first get the list of confident NIMs from [Sankararaman et al 2014](https://www.nature.com/articles/nature12961) in 1KG EUR populations. Then we expand it to include all the SNPs in strong LD (r^2 >= 0.99) with any confident NIMs as expanded NIMs with `plink --bfile <qced>  --show-tags <confident NIMs>  --tag-r2 0.99  --tag-kb 200 --out <expanded NIMs>`

 ```diff
 - note to self: 
 all results for tag SNPs and tag SNPs annotation are based on 'ND.original.id' so everything downstream from here does not need any update.
 ```

Calculate if the tag SNPs are more likely to match with Altai Neanderthal reference with `crfTags.m` where the CRF summaries need to be downloaded from %%% and the summaries used are in `summaries.release/EUR.hapmap/summaries/` 

```diff
- note to self: 
`crfTags.m` in  `/u/home/s/sriram/group/sriram/projects/ukbio/april/nimHeritability/simulation/adaptedScripts/match
```

### MAF
We use 5 MAF bins, which split all QC-ed SNPs from low to high MAF into equal sized bins
### LD
We use 5 LD bins, which split all QC-ed SNPs from low to high LD-score into equal sized bins.

### Consruct non-overlapping annotation

In the output file name annotation for ancestry is indicated as '.anc', maf with '.maf', and ld-score with '.ld'. 

The annotation is non-overlapping, meaning that SNPs are partitioned into combinations of annotation, hence each SNP belongs to one and only one annotation. 

For example, in a file named with 'nol.anc.maf.annot', there are 10 annotations total (2 ancestry * 5 maf), with the Neanderthal ancestry into 5 maf bins, and modern human ancestry into 5 maf bins. 

Annotation  with confident NIM: `writeAnnotOri.m`, output 'allOri.nol.anc.annot', 'allOri.nol.anc.maf.annot', 'allOri.nol.anc.ld.annot'

Annotation with tag SNPs: `writeAnnotTags.m`, output 'tags.99.nol.anc.annot', 'tags.99.nol.anc.maf.annot', 'tags.99.nol.anc.ld.annot'

Ancestry-MAF-LD based Annotation with tag SNPs: `writeAnnotMAFLD.m`, output 'tags.99.nol.anc.maf.ld.annot', and 'indNearEmptyBins.tags.nol.maf.ld.txt'. 

*An exception is '.anc.maf.ld' annotation which should have 2*5*5 = 50 non-overlapping bins, but because there are few high MAF Neanderthal SNPs and low LD-score Neanderthal SNPs, some bins are empty or near empty, hence we remove the bins with smaller than 30 SNPs in them. This particular annotation is only constructed for Tag SNPs because the number of confident NIMs is too small to have so many annotations.*
```diff
- note to self: 
`/u/home/s/sriram/group/sriram/projects/ukbio/april/nimHeritability/simulation/adaptedScripts/match
```
## Estimating heritability with RHE-mc
### Information about RHE-mc
Information about RHE-mc can be found at [RHE-mc GitHub](https://github.com/sriramlab/RHE-mc)
### Simulated data
The gcta64 simulated phenotype does not have a header, e.g.:

        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 
        1000075 1000075 60.4556
This file can be directly used with plink for GWAS but RHE-mc takes phenotype file which starts with a header, hence require adding a header line with

        `sed  -i '1 i\FID IID pheno' *.phen`
Then we can get the correct format for RHE-mc, e.g.:

        FID IID pheno
        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 
        1000075 1000075 60.4556
Run RHE-mc with the supply of genotype, phenotype, and annotation files for whole-genome simulated data with

        RHEmc -g <genotype file> (e.g., qced)  -p <phenotype file> -annot <annotation file> -k 10 -jn 100  -o <output file>

### UKBB data       
Run RHE-mc with the supply of genotype, phenotype, coavariate and anonotation files for UKBB data.

        RHEmc -g <genotype file> (e.g., qced)  -p <phenotype file> -c <covar file>  -annot <annotation file> -k 10 -jn 100  -o <output file>
   
### H2 Partitioning
From the output of RHE-mc, we extract the heritabiliity and standard error of heritability for each annotation in the non-overlapping setting to get the partitioned heritability estimates for SNPs in Neanderthal ancestry.

```diff
- note to self: 
currently in personal computer **and need to reanalyze these with the allOri.annot instead of the ori.annot
```
### UKBB META-analysis

```diff
- note to self: 
code on computer `metaAnalysis.m`
```

## Fine mapping

### Information about SuSiE 
Information about Susie can be found at [SuSiE GitHub](https://stephenslab.github.io/susie-paper/index.html). It is an R package for fine mapping.
        
        
### Benchmark SuSiE with simulated data

In this part, we perform from three different sets of tests: GWAS significant SNPs, GWAS significant NIMs, sliding window across the genome.
#### GWAS significant SNPs
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
4. Analyze all Susie output with `getSusieStats.m`

#### GWAS significant NIMs
1. GWAS with plink, then extrat significant NIMs with P < 10^(-10) with `getSignificantSNPs.m`
2. For each GWAS signficant NIM
        2.1 Output the 200 kb region surrounding this NIM
        2.2 Run SuSiE on the 200kb region surrounding the NIM with script `runSusie.R` 
3. Analyze all Susie output with `getSusieStats.m`

#### Sliding window analysis

1. Divide genomes into 200kb windows with 100kb step size.
2. For each window:
        2.1 Output the 200 kb region 
        2.2 Run SuSiE on the 200kb region with script `runSusie.R` 
3. Analyze all Susie output with `getSusieStats.m`

```diff
- note to self: 
/u/home/s/sriram/group/sriram/projects/ukbio/april/nimHeritability/simulation/adaptedScripts/fineMapping
```
### SuSiE to UKBB phenotypes

```diff
- note to self:
Collection of phenotypes: https://docs.google.com/spreadsheets/d/1L7XNxSm1M8rGbJy0G8TydsDhiIgzImhs9LUOZxgN8XM/edit?usp=sharing
```
