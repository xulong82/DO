# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: show the results
# Rev: April 6, 2014
#---------------------------------------------------------------------------------------------------
library(ggplot2)
library(parallel)
#---------------------------------------------------------------------------------------------------
# load("~/Dropbox/DO/R/myanalysis1.rdt")
# source("~/Dropbox/DO/R/myinfo.R")  # chrlen, chrmid, ens.ucsc, snp.inf
load("~/Dropbox/DO/R/myinfo.rdt")
# load("~/HPC/doqtl.rdt")  # DOQTL
# load("~/HPC/doqtl.perm1.rdt")  # DOQTL
# load("~/HPC/doqtl.perm2.rdt")  # DOQTL
# load("~/HPC/doqtl.perm3.rdt")  # DOQTL
# save(doqtl.perm, file = "~/Dropbox/DO/R/doqtl.perm.rdt")
load("~/Dropbox/DO/R/doqtl.perm.rdt")  # DOQTL
#---------------------------------------------------------------------------------------------------
lod.thr <- quantile(doqtl.perm, 0.95)

#---------------------------------------------------------------------------------------------------
# gene.ind <- which(rownames(pheno) %in% rownames(ens.ucsc))
# gene.ids <- rownames(pheno)[which(rownames(pheno) %in% rownames(ens.ucsc))]
# 
# doqtl.list <- doqtl.list[gene.ind]
# n.gene <- length(gene.ind)
# pos.gene <- ens.ucsc[gene.ids, ]$pos
# 
# snp.a <- doqtl.list[[1]]$lod$A$SNP_ID  # SNPs in autosome
# snp.x <- doqtl.list[[1]]$lod$X$SNP_ID  # SNPs in X-chromosome
# # snp.ids <- c(snp.a, snp.x) 
# snp.ids <- snp.a
# pos.snp <- snp.inf[snp.ids, ]$pos
# null <- which(is.na(pos.snp))
# snp.ids <- snp.ids[-null]
# n.geno <- length(snp.ids)
# pos.snp <- pos.snp[-null]
# save(gene.ids, pos.gene, n.gene, snp.ids, pos.snp, n.geno, file = "~/Dropbox/DO/R/myinfo2.rdt")
load("~/Dropbox/DO/R/myinfo2.rdt")

#---------------------------------------------------------------------------------------------------
# doqtl.lod <- matrix(nrow = n.gene, ncol = n.geno, dimnames = list(gene.ids, snp.ids))
# coef1 <- colnames(doqtl.list[[1]]$coef$A)
# doqtl.coef <- array(dim = c(n.gene, n.geno, length(coef1)), dimnames = list(gene.ids, snp.ids, coef1))
# for (i in 1:n.gene) {
#   if (i %% 100 == 0) print(i)
# # doqtl.lod[i, ] <- doqtl.list[[i]]$lod$A$lod[-null]
#   doqtl.coef[i, , ] <- doqtl.list[[i]]$coef$A[-null, ]
# }
load("~/HPC/doqtl.lod.rdt")
load("~/HPC/doqtl.coef.rdt")

hist.dt1 <- data.frame(lod = as.vector(doqtl.lod))
pdf("~/Dropbox/DO/Figures/hist1.pdf")
hist(doqtl.lod, freq = F, main = "", xlab = "LOD", 
     font = 2, font.lab = 2, nclass = 1500, xlim = c(0, 10), lwd = 2)
abline(v = 7.31, col = "black", lwd = 2)
dev.off()

#---------------------------------------------------------------------------------------------------
# qtl.lod1 <- data.frame(row.names = snp.ids, cis = rep(F, n.geno), trans = rep(F, n.geno)) 
# qtl.lod2 <- 0
# 
# for (i in 1:n.geno) {  # Loop over SNPs
#   if (i %% 100 == 0) print(i)
#   lod <- doqtl.lod[, i]
#   pos <- pos.snp[i]  # SNP position
#   for (j in 1:n.gene) {
#     if (lod[j] > lod.thr) {
#       qtl.lod2 <- c(qtl.lod2, c(i, j, lod[j]))
#       if (abs(pos - pos.gene[j]) < 30) {
#         qtl.lod1[i, 1] = T
#       } else {
#         qtl.lod1[i, 2] = T
#       }
#     }
#   }
# }
load("~/HPC/qtl.lod1.rdt")
load("~/HPC/qtl.lod2.rdt")

#---------------------------------------------------------------------------------------------------
qtl.lod1.1 <- qtl.lod1[qtl.lod1$cis == T & qtl.lod1$trans == F, ] # cis
qtl.lod1.2 <- qtl.lod1[qtl.lod1$cis == F & qtl.lod1$trans == T, ] # trans
qtl.lod1.3 <- qtl.lod1[qtl.lod1$cis == T & qtl.lod1$trans == T, ] # cis & trans
cis <- which(snp.ids %in% rownames(qtl.lod1.1))
trans <- which(snp.ids %in% rownames(qtl.lod1.2))
cis.trans <- which(snp.ids %in% rownames(qtl.lod1.3))

qtl.lod2 <- matrix(qtl.lod2[-1], ncol = 3, byrow = T)
colnames(qtl.lod2) <- c("snp_id", "gene_id", "lod")
qtl.lod2.1 <- qtl.lod2[qtl.lod2[, 1] %in% cis, ]
qtl.lod2.2 <- qtl.lod2[qtl.lod2[, 1] %in% trans, ]
qtl.lod2.3 <- qtl.lod2[qtl.lod2[, 1] %in% cis.trans, ]

str(unique(qtl.lod2[, 1]))
str(unique(qtl.lod2.1[, 1]))
str(unique(qtl.lod2.2[, 1]))
str(unique(qtl.lod2.3[, 1]))

mysnps <- rownames(qtl.lod1.3)
# mysnps <- mysnps[snp.inf[mysnps, ]$Chr != 20]  # Sex chromosome
myqtls <- doqtl.lod[, mysnps]

# ens2sym <- read.delim("~/Dropbox/Genome/ensembl2symbol.map", header = F)
# ens2sym <- data.frame(row.names = ens2sym$V1, Symbol = ens2sym$V2)
# mygenes.sym <- ens2sym[mygenes, ]
# tfs <- read.delim("~/Dropbox/Genome/TFdb.Riken.txt", header = F)  # 1675 TFs from Riken TFdb
# mytfs <- intersect(tfs$V1, gene.sel.sym)  # 900 expressed TFs

map.dt1 <- data.frame(pos.snp = pos.snp[qtl.lod2.3[, 1]], pos.gene = pos.gene[qtl.lod2.3[, 2]], lod = qtl.lod2.3[, 3])
pdf("~/Dropbox/DO/Figures/map1.pdf")
ggplot(map.dt1, aes(x = pos.snp, y = pos.gene)) + geom_point(aes(size = lod)) + 
  theme_bw() + 
  xlab("SNP") + ylab("Gene") +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
  scale_y_continuous(breaks = chrmid, labels = names(chrlen)) +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold")) +
#       axis.title.x = element_text(vjust = -0.5)) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

# myMap <- function (map.dt) {
#   plot(map.dt[, 1:2], pch = 20, cex = .2, xlim = c(0, max(chrlen)), ylim = c(0, max(chrlen)), xaxt = "no", yaxt = "no")
#   abline(h = chrlen, col = "grey80")
#   abline(v = chrlen, col = "grey80")
#   axis(1, at = chrmid, names(chrlen))
#   axis(2, at = chrmid, names(chrlen))  
# }
