# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: DO:DOQTL/myQtl1
# Rev: April 9, 2014
#---------------------------------------------------------------------------------------------------
library(DOQTL)
library(parallel)

setwd("~/Dropbox/DO")

load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))  # MUGA annotation
load("Nara/CheslerDO.founder.probs.mapping.Rdata")  # founder.probs | DOQTL:calc.genoprob
# kinship.nara <- read.delim("Nara/CheslerDO.K.mapping.txt")  # Nara's Kinship matrix
pheno <- read.delim("Nara/CheslerDOexp-se-merged-rzqnorm.txt")
colnames(pheno) <- gsub("[MF]", "", colnames(pheno))
batch <- read.delim("Nara/CheslerDO.batch.mapping.txt")

founder.probs <- founder.probs[sort(rownames(founder.probs)), , ]
pheno <- pheno[, sort(colnames(pheno))]
kinship <- kinship.probs(founder.probs)  # DOQTL:kinship.probs

covar <- data.frame(row.names = batch$mID, sex = batch$sex, run = batch$runD1)
covar <- covar[sort(rownames(covar)), ]
covar <- data.frame(row.names = rownames(covar), covar$sex)

n.sample <- dim(pheno)[2]
n.pheno <- dim(pheno)[1]
n.geno <- dim(founder.probs)[3]

pheno.list <- list()
for (i in 1:nrow(pheno)) pheno.list[[i]] <- t(pheno[i, ])
founder.probs.list <- list()
for (i in 1:ncol(founder.probs[1, , ])) founder.probs.list[[i]] <- founder.probs[, , i]

save.image("~/Dropbox/DO/R/myanalysis1.rdt")
load("~/Dropbox/DO/R/myanalysis1.rdt")

#-- Check the data ----------------------------------------------------------------------------------
if ("FALSE" %in% (dimnames(table(rownames(covar) == rownames(founder.probs[, , 1]))))) {
  stop("Samples do not match!")
}
if ("FALSE" %in% (dimnames(table(colnames(pheno) == rownames(founder.probs[, , 1]))))) {
  stop("Samples do not match!")
}

#-- DOQTL ------------------------------------------------------------------------------------------
mc.scanone <- function (x) {  # single run
  y <- scanone(pheno = x, probs = founder.probs, K = kinship, addcovar = covar, snps = muga_snps)
# return(rbind(y$lod$A, y$lod$X))
  return(y)
}
mc.scanone.perm <- function (x) {  # permutation
  y1 <- permutations.qtl.LRS(pheno = x, probs = founder.probs, addcovar = covar, 
                             snps = muga_snps, nperm = 100)
  return(y1 / (2 * log(10)))  # Lod
}

#-- Parallel computing ----------------------------------------------------------------------------
time.start <- strptime(date(), "%a %b %d %H:%M:%S %Y")
# doqtl.list <- mclapply(pheno.list, mc.scanone, mc.cores = 5)
doqtl.perm.list <- mclapply(pheno.list[1:5], mc.scanone.perm, mc.cores = 5)
time.end <- strptime(date(), "%a %b %d %H:%M:%S %Y")
difftime(time.end, time.start, units = "auto")

# #-- MYQTL:lm ---------------------------------------------------------------------------------------
# myQtl1 <- function (x) {  # R:lm() without kinship adjustment
#   # print(colnames(pheno1))  # Progress
#   data0 <- data.frame(pheno = x[, 1], covar = covar[, 1])
#   g0 <- lm(pheno ~ ., data = data0)
#   y <- rep(0, n.geno)
#   for (j in 1:n.geno) {
#     data1 <- data.frame(pheno = x[, 1], covar = covar[, 1], founder.probs[, , j])
#     g1 <- lm(pheno ~ ., data = data1)
#     lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
#     lod <- lrt / (2 * log(10))  # Lod
#     y[j] <- lod
#   }
#   return(y)
# }
# myQtl1.perm <- function (x) {  # R:lm() without kinship adjustment
#   # print(colnames(pheno1))  # Progress
#   y = 1
#   for (i in 1:n.perm) {  # permutation times
#     x1 <- matrix(sample(x))
#     dimnames(x1) <- dimnames(x)
#     data0 <- data.frame(pheno = x1[, 1], covar = covar[, 1])
#     g0 <- lm(pheno ~ ., data = data0)
#     #   print(i)
#     y1 <- rep(0, n.geno)
#     for (j in 1:n.geno) {
#       data1 <- data.frame(pheno = x1[, 1], covar = covar[, 1], founder.probs[, , j])
#       g1 <- lm(pheno ~ ., data = data1)
#       lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
#       lod <- lrt / (2 * log(10))  # Lod
#       y1[j] <- lod
#     }
#     y <- c(y, y1)
#   }
#   return(y[-1])
# }
# #-- QTLrel:estVC and lm ----------------------------------------------------------------------------
# myQtl2 <- function (pheno1) {  # R:lm() with kinship adjustment
#   # pheno1 <- t(pheno[i, ])
#   print(colnames(pheno1))  # Progress
#   
#   VC <- estVC(pheno1[, 1], covar, 
#               v = list(AA = kinship, DD = NULL, HH = NULL, AD = NULL, MH = NULL, EE = diag(nrow(kinship))))
#   V <- VC$par["AA"] * kinship + VC$par["EE"] * diag(nrow(kinship))
#   eig <- eigen(V)
#   A <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
#   
#   pheno1.t <- solve(A) %*% pheno1
#   covar.t <- solve(A) %*% as.matrix(covar)
#   
#   data0 <- data.frame(pheno = pheno1.t[, 1], covar = covar.t[, 1])
#   g0 <- lm(pheno ~ ., data = data0)
#   
#   myqtl2 <- data.frame(lrt = rep(0, n.geno), lod = rep(0, n.geno))
#   for (j in 1:n.geno) {
#     data1 <- data.frame(pheno = pheno1.t[, 1], covar = covar.t[, 1], founder.probs[, , j])
#     g1 <- lm(pheno ~ ., data = data1)
#     
#     lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
#     lod <- lrt / (2 * log(10))  # Lod
#     myqtl2[j, ] <- c(lrt, lod)
#   }
#   return(myqtl2)
# }
