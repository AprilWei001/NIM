# NIM

## Table of contents
* [UKBB SNP QC](#UKBB-SNP-QC)
* [Whole-genome simulation](#Whole-genome-simulation)
* [SNP annotation](#SNP-annotation)
* [Estimating heritability with RHE-mc](#Estimating-heritability-with-RHE-mc)
* [Fine mapping](#Fine-mapping)


## UKBB SNP QC 
SNP QC for UKBB whole-genome imputed data with plink flags 
        
       plink  --bfile <ukbb imputed genotype data> --exclude <SNPs in the mhc region txt file> --maf 0.001 --geno 0.01 --hwe 1e-7 --keep <ukbb wb unrelated individuals fam file> --make-bed --out <output file name (e.g., qced)>


```diff
-note to self: 
find this information on hoffman /u/project/sgss/UKBB/data/imp/Qced/All/README, although the QC-ed data used are in /u/project/sgss/UKBB/data/imp/Qced/maf_0.001 which persumably only changes the flag --maf 0.01 into --maf 0.001
```

## Whole-genome simulation
Script: `simAnyArchitecture.sh`
* This simulation depends on software: plink, gcta64, cripts: `shuffle.py`
* input files: <genotype input file name (e.g. qced)>  <.frq> which has the in-sample MAF of SNPs used, and <.ld> in-sample ldscore: comptued with 

      `gcta64 --bfile <genotype input file name (e.g. qced)> --ld-score --ld-wind 10000 --out <output filename>` 

Run with:

    `./simAllArchitecture.sh <type of architecture> <seed>`
for example 

    `./simAllArchitecture.sh POLY 1` 
will give the simulation for polygenic architecture with seed 1. 

The types of architecture can be:
* POLY -- polygenic architecture 
* LOW -- low LD SNPs (defined as ldsc <=10) contributes to 9000 causal variants, and high LD SNPs contributes 1000 causal variants, 
* HIGH -- high LD (defined as ldsc > 10) contributes to 9000 causal variants, and low LD SNPs contributes 1000 causal variants,
* RARE -- rare SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and common SNPs contributes 1000 causal variants, 
* COMMON -- common SNPs (defined as maf <= 0.05) contributes to 9000 causal variants, and rare SNPs contributes 1000 causal variants,


```diff
-note to self: 
script and small sized dependencies moved to a folder on tabla: /home/aprilwei/projects/nimHeretability/github -- original simulations are outputs from other scripts but merged here for the purpose of sharing, tested runnable on tabla.
```

## SNP annotation
### Ancestry 
We use two annotation of ancestry: neanderthal ancestry and modern human ancestry. Neanderthal ancestry is defined either with confident NIM or with tag SNPs, and modern human ancestry is therefore the SNPs that are not confident NIMs, or the SNPs that are not tag SNPs.

* confident NIM

Extract the coordinate id from `nim.confident.bim` with 

        `cat /u/home/s/sriram/group/sriram/projects/ukbio/data/geno/imp/filtered/nim.confident.bim |awk '{print $2}' > ND.original.id`

```diff
- note to self: 
originally also did `intersect -k1 0 -k2 0 -f1 ND.original.id -f2 UKBB.maf_daf_ld > ND.id` to get an interesect for SNP-by-SNP matching, and these are used in `writeAnnot.m` but then was updated with `writeAnnotOri.m` on Jan 14 and used in downstream anayses, so everything is okay!
```

* tag SNPs

Identify UKBB SNPs in high LD with confident NIMs with `plink --bfile <qced>  --show-tags ND.original.id  --tag-r2 0.99  --tag-kb 200 --out nimLD.99`

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

### generate non-overlapping annotation

In the output file name annotation for ancestry is indicated as '.anc', maf with '.maf', and ld-score with '.ld'. 

The annotation is non-overlapping, meaning that SNPs are partitioned into combinations of annotation, hence each SNP belongs to one and only one annotation. 

For example, in a file named with 'nol.anc.maf.annot', there are 10 annotations total (2 ancestry * 5 maf), with the Neanderthal ancestry into 5 maf bins, and modern human ancestry into 5 maf bins. *An exception is '.anc.maf.ld' annotation which should have 2*5*5 = 50 non-overlapping bins, but because there are few high MAF Neanderthal SNPs and low LD-score Neanderthal SNPs, some bins are empty or near empty, hence we remove the bins with smaller than 30 SNPs in them.

Annotation  with confident NIM: `writeAnnotOri.m`, output 'allOri.nol.anc.annot', 'allOri.nol.anc.maf.annot', 'allOri.nol.anc.ld.annot'

Annotation with tag SNPs: `writeAnnotTags.m`, output 'tags.99.nol.anc.annot', 'tags.99.nol.anc.maf.annot', 'tags.99.nol.anc.ld.annot'

Ancestry-MAF-LD based Annotation with tag SNPs: `writeAnnotMAFLD.m`, output 'tags.99.nol.anc.maf.ld.annot', and 'indNearEmptyBins.tags.nol.maf.ld.txt'. 

```diff
- note to self: 
`/u/home/s/sriram/group/sriram/projects/ukbio/april/nimHeritability/simulation/adaptedScripts/match
```
## Estimating heritability with RHE-mc
### Simulated data
The simulated phenotype with gcta64 does not have a header, for example:

        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 
        1000075 1000075 60.4556
This file can be directly used with plink for GWAS but RHE-mc takes phenotype file which starts with a header, hence require adding a header line with

        `sed  -i '1 i\FID IID pheno' *.phen`
Then we can get

        FID IID pheno
        1000026 1000026 -179.62 
        1000058 1000058 -116.748 
        1000060 1000060 123.494 
        1000075 1000075 60.4556
Run RHE-mc with the supply of genotype, phenotype, and annotation files for whole-genome simulated data with

        `RHEmc -g <genotype file> (e.g., qced)  -p <phenotype file> -annot <annotation file> -k 10 -jn 100  -o <output file>`
        
### UKBB data       
Run RHE-mc with the supply of genotype, phenotype, coavariate and anonotation files for UKBB data.

        `RHEmc -g <genotype file> (e.g., qced)  -p <phenotype file> -c <covar file>  -annot <annotation file> -k 10 -jn 100  -o <output file>`
   
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
### GWAS
Use plink with flag '--linear standard-beta'
        
        `plink --silent --bfile <qced>--pheno $OUT.phen --linear standard-beta --out <output file>`
### Susie 
Information about Susie can be found at (https://stephenslab.github.io/susie-paper/index.html)

After installing Susie in R, one can perform a single locus fine mapping which outputs confidence sets with script `runSusie.R`
Run with:

        R --slave --args <genotype text file > <phenotype file> < runSusie.R >
        
This rscript takes the genotype text file input from the tested region. Such genotype file can be generated from 'bfile' with:

        plink --bfile <qced> --snps <range of the snps (e.g. 1:4669624-1:4869948)> --recode A --out <output genotype text file>
        
### Benchmark Susie with simulated data

```diff
- note to self: 
/u/home/s/sriram/group/sriram/projects/ukbio/april/nimHeritability/simulation/adaptedScripts/fineMapping
```
### Apply to UKBB phenotypes
