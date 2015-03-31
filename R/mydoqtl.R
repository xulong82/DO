# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Diversity of outbred (DO) analysis
# Rev: March 19, 2014

library(DOQTL)
library(parallel)

setwd("~/Dropbox/DO")

#---------------------------------------------------------------------------
load("./MUGAExampleData/data/FinalReport1.rda")
load("./MUGAExampleData/data/FinalReport2.rda")
load("./MUGAExampleData/data/Samples1.rda")
load("./MUGAExampleData/data/Samples2.rda")
load("./MUGAExampleData/data/pheno.rda")
load("./MUGAExampleData/data/geno.rda")
load("./MUGAExampleData/data/x.rda")
load("./MUGAExampleData/data/y.rda")
load("./MUGAExampleData/data/model.probs.rda")
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))  # MUGA annotation

K <- kinship.probs(model.probs)  # Calculate the kinship matrix using the founder contributions
pheno <- pheno[rownames(pheno) %in% rownames(K), ]
covar <- data.frame(row.names = rownames(pheno),
                    sex = as.numeric(pheno$Sex == "M"), diet = as.numeric(pheno$Diet == "hf"))

pheno$HDW2[is.na(pheno$HDW2)] <- 0

qtl <- scanone(pheno = pheno, pheno.col = "HDW2", probs = model.probs, 
               K = K, addcovar = covar, snps = muga_snps)
perms <- scanone.perm(pheno = pheno, pheno.col = "HDW2", probs = model.probs, 
                      K = K, addcovar = covar, snps = muga_snps, nperm = 10)

interval <- bayesint(qtl, chr = 9)

ma <- assoc.map(pheno = pheno, pheno.col = "HDW2", probs= model.probs,
                K = K, addcovar = covar, snps = muga_snps, chr = interval[1, 2],
                start = interval[1, 3], end = interval[3, 3])

save.image("~/Dropbox/DO/myanalysis1.20140320.RData")
load("~/Dropbox/DO/myanalysis1.20140320.RData")
