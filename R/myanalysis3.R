# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: DO:PCA
# Rev: March 19, 2014
#---------------------------------------------------------------------------------------------------
load("~/Dropbox/DO/R/allDO.genotypes.Rdata")
geno.ac <- geno[grep("chesler", rownames(geno)), ]  # genotype in allele call
rownames(geno.ac) <- gsub("chesler", "", rownames(geno.ac))
rm(geno)

myRef <- function (x) {
  y1 <- as.matrix(table(x))
  y2 <- intersect(rownames(y1), c("A", "T", "C", "G"))
  if (length(y2) == 0) {
    y = "U"
  } else {
    y <- names(sort(y1[y2, ], decreasing = T))[1]
  }
  return(y)
}

refs <- apply(geno.ac, 2, myRef)

geno.bi <- matrix(nrow = nrow(geno.ac), ncol = ncol(geno.ac))
dimnames(geno.bi) <- dimnames(geno.ac)
for (i in 1:nrow(geno.ac)) {
  for (j in 1:ncol(geno.ac)) {
    if (geno.ac[i, j] == refs[j]) {
      geno.bi[i, j] <- 0
    } else if (geno.ac[i, j] == "H") {
      geno.bi[i, j] <- 1
    } else if (geno.ac[i, j] == "N") {
      geno.bi[i, j] <- geno.bi[i, j]
    } else {
      geno.bi[i, j] <- 2
    }
  }
}
#---------------------------------------------------------------------------------------------------
load("~/Dropbox/DO/R/myanalysis1.rdt")
geno <- geno.bi[colnames(pheno), ]

load("~/Dropbox/DO/R/mysnps.rdt")
sanger <- mysnps[intersect(rownames(mysnps), snp.ids), ]
probs <- founder.probs[, , rownames(sanger)]
n.geno2 <- dim(probs)[3]

geno.ac <- matrix(0, nrow = n.sample, ncol = n.geno2)
for (i in 1:n.sample) {
  for (j in 1:n.geno2) {
    geno.ac[i, j] <- sanger[j, 1:8] %*% probs[i, 1:8, j]
  }
}
rownames(geno.ac) <- rownames(probs[, , 1])
colnames(geno.ac) <- colnames(probs[1, , ])

pheno1 <- pheno[gene.sel1, ]
#-- PCA and lm -------------------------------------------------------------------------------------
svd1 <- prcomp(geno.ac)  # PC in terms of genotype
barplot(svd1$sdev, xlab = "PC", ylab = "Standard Deviation")
abline(h = 0)

pcs <- svd1$x[, 1:5]  # predict(svd1)

myQtl3 <- function (x) {  # R:lm() with PCA correction
  # pheno1 <- t(pheno[i, ])
  print(rownames(x))  # Progress
  
  data <- data.frame(pheno = t(x)[, 1], pcs)
  pheno1.t <- lm(as.formula(paste("pheno", "~",
                 paste(colnames(data)[c(2:ncol(data))], collapse = "+"), sep = "")), data)$res
  
  data0 <- data.frame(pheno = pheno1.t, covar = covar1)
  g0 <- lm(pheno ~ ., data = data0)
  
  y <- rep(0, n.geno2)
  for (j in 1:n.geno2) {
    data1 <- data.frame(pheno = pheno1.t, covar = covar1, probs[, , j])
    g1 <- lm(pheno ~ ., data = data1)
    
    lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
    lod <- lrt / (2 * log(10))  # Lod
    y[j] <- lod
  }
  return(y)
}
#---------------------------------------------------------------------------------------------------
pheno.list <- list()  # Test purpose
for (i in 1:10) pheno.list[[i]] <- t(pheno[i, ])
gene.tst <- pheno1[6, ]
lod.tst <- myQtl3(gene.tst)

pos.snp2 <- pos.snp[which(snp.ids %in% rownames(sanger))]
par(mfrow = c(3, 1))
plot(pos.snp, qtl.sel[8, ], type = "l", main = "DOQTL")
plot(pos.snp, myqtl.sel.lod[6, ], type = "l", main = "MYQTL")
plot(pos.snp2, lod.tst, type = "l", main = "PCA")

themax.list <- t(pheno["ENSMUSG00000020268", ])
themax.qtl = myQtl3(themax.list)

