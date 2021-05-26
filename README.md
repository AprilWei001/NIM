# NIM
## UKBB SNP QC 
QC with plink flags `--exclude <mhc region> --maf 0.001 --geno 0.01 --hwe 1e-7 --keep <ukbb wb unrelated individuals> --make-bed --out <qced data>`

## Whole genome simulation
for any genetic architecture on tabla: /home/aprilwei/projects/nimHeretability/wgSim/simAllArchitecture.sh
This file takes UKBB in-sample allele frequency (for simulating MAF based architecture) and in-sample ldscore (for simulating LD based architecture) as input.
