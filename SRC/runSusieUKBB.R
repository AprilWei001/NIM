library(susieR)

arg = commandArgs(trailingOnly=T)
filename = arg[1]
phenoname = arg[2]

A <- read.table(filename, header = TRUE, sep = " ", dec = ".")

pheno <- read.table(phenoname, header = TRUE, sep = " ", dec = ".")

pheno = pheno[[3]]#third col is the phenotypic values, where -9 is the same as na

pheno[pheno == -9]<- NA

b <- dim(A)

A = A[, 7:b[2]]

normalize = function(x){
if(length(which(is.na(x)))==0) (x-mean(x))/sd(x) else
(x-mean(x,na.rm=T))/sd(x,na.rm=T)
}

A = apply(A, 2, normalize)

A[is.na(A)] = 0
 
A = A[is.na(pheno) == FALSE,]
pheno = pheno[is.na(pheno) == FALSE]

fitted <- susie(A, pheno,
                L = 10,
                estimate_residual_variance = TRUE, 
                estimate_prior_variance = FALSE,
                scaled_prior_variance = 0.1,
                verbose = TRUE)
index <- fitted$sets$cs_index

snpID <- colnames(A)
 
print(fitted$sets)

for (ind in index)
{
	print(snpID[fitted$sets$cs[[ind]]])
}

