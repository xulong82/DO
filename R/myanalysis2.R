# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: DO:cape/capeDO
# Rev: March 19, 2014

library(cape)
library(parallel)

setwd("~/Dropbox/DO")
load("~/Dropbox/DO/R/myanalysis1.rdt")
load("~/Dropbox/DO/R/doqtl1.rdt")
source("~/Dropbox/DO/R/myannots.R")  # chrlen, chrmid, ens.ucsc, snp.inf

library(capeDO)
# library(capeAddons)

cross <- list()
cross$geno <- founder.probs[, , snp.sel]
cross$pheno <- t(pheno)[, gene.sel1]

cross$marker.names <- snp.sel
cross$chromosome <- snp.inf[snp.sel, ]$Chr
cross$location <- snp.inf[snp.sel, ]$Mb_NCBI38
cross$snp.pos <- snp.inf[snp.sel, ]$pos 
cross$non.allelic.covar <- covar

cross <- get.eigentraits(cross, scale.pheno = FALSE, normalize.pheno = FALSE)
cross <- select.eigentraits(cross, traits.which = c(1:10))
plotSVD(cross, orientation = "vertical")

names(dimnames(cross$geno)) <- c("mouse", "allele", "locus")
cross <- singlescan(cross, n.perm = 10, ref.allele = "B", covar = "covar.sex", 
                    alpha.for.pairs = 0.05, alpha.for.covar = 0.05, 
                    scan.what = "eigentrait", auto.covar.selection = FALSE, verbose = TRUE)

covar1 <- covar[, 1]
n.snp1 <- length(snp.sel)
cape.lod <- matrix(nrow = 10, ncol = n.snp1)
for (i in 1:10) {
  print(colnames(cross$ET)[i])  # Progress
  x <- cross$ET[, i]
  data0 <- data.frame(pheno = x, covar = covar1)
  g0 <- lm(pheno ~ ., data = data0)
  lod <- rep(0, n.snp1)
  for (j in 1:n.snp1) {
    data1 <- data.frame(pheno = x, covar = covar1, cross$geno[, , j])
    g1 <- lm(pheno ~ ., data = data1)
    lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
    lod[j] <- lrt / (2 * log(10))  # Lod
  }
  cape.lod[i, ] <- lod
}

pos.snp1 <- pos.snp[which(snp.ids %in% snp.sel)]
plot(pos.snp1, cape.lod[10, ], type = "p", pch = 20, cex = .2)

line.dt2 <- data.frame(pos = 0, lod = 0, gene = 0)
for (i in 1:10) {
  tmp <- data.frame(pos = pos.snp1, lod = cape.lod[i, ], gene = rep(paste("ET", i, sep = ""), n.snp1))
  line.dt2 <- rbind(line.dt2, tmp)
}
line.dt2 <- line.dt2[-1, ]

ggplot(line.dt2, aes(x = pos, y = lod, group = gene, colour = factor(gene))) +
  geom_point() + 
  geom_abline(intercept = 8, slope = 0) +
  theme_bw()

cape.lod.perm <- list()
for (i in 1:10) {
  print(colnames(cross$ET)[i])  # Progress
  x <- cross$ET[, i]
  y = 0
  for (k in 1:1000) {  # permutation times
    print(k)
    x1 <- sample(x)
    data0 <- data.frame(pheno = x1, covar = covar1)
    g0 <- lm(pheno ~ ., data = data0)
    y1 <- rep(0, 2072)
    for (j in 1:2072) {
      data1 <- data.frame(pheno = x1, covar = covar1, cross$geno[, , j])
      g1 <- lm(pheno ~ ., data = data1)
      lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
      lod <- lrt / (2 * log(10))  # Lod
      y1[j] <- lod
    }
    y <- c(y, y1)
  }
  cape.lod.perm[[i]] <- y[-1]
}

saveRDS(cross, "crossDO.RData")
