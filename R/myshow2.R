library(gplots)

load("~/HPC/myqtl1.rdt")  # DOQTL

myqtl1.list <- myqtl1.list[gene.ind]

myqtl1.lod <- matrix(nrow = n.gene, ncol = n.geno)
for (i in 1:n.gene) {
  if (i %% 1000 == 0) print(i)
  myqtl1.lod[i, ] <- myqtl1.list[[i]]
}
myqtl1.lod <- myqtl1.lod[, -null1]

myqtl1 <- data.frame(row.names = gene.ids, cis = rep(F, n.gene), trans = rep(F, n.gene)) 
myqtl2 <- 0

for (i in 1:n.gene) {
  if (i %% 100 == 0) print(i)
  lod <- myqtl1.lod[i, ]
  pos <- pos.gene[i]
  for(j in 1:n.geno1) {
    if (lod[j] > lod.thr) {
      myqtl2 <- c(myqtl2, c(i, j, lod[j]))
      if (abs(pos - pos.snp[j]) < 30) {
        myqtl1[i, 1] = T
      } else {
        myqtl1[i, 2] = T
      }
    }
  }
}

myqtl1.1 <- myqtl1[myqtl1$cis == T & myqtl1$trans == T, ]  # cis and trans
myqtl1.2 <- myqtl1[myqtl1$cis == T & myqtl1$trans == F, ]  # cis
myqtl1.3 <- myqtl1[myqtl1$cis == F & myqtl1$trans == T, ]  # trans
mycis.trans <- which(gene.ids %in% rownames(myqtl1.1))
mycis <- which(gene.ids %in% rownames(myqtl1.2))
mytrans <- which(gene.ids %in% rownames(myqtl1.3))

myqtl2 <- matrix(myqtl2[-1], ncol = 3, byrow = T)
colnames(myqtl2) <- c("gene_id", "snp_id", "lod")
myqtl2.1 <- myqtl2[myqtl2[, 1] %in% mycis.trans, ]
myqtl2.2 <- myqtl2[myqtl2[, 1] %in% mycis, ]
myqtl2.3 <- myqtl2[myqtl2[, 1] %in% mytrans, ]

str(unique(myqtl2[, 1]))
str(unique(myqtl2.1[, 1]))
str(unique(myqtl2.2[, 1]))
str(unique(myqtl2.3[, 1]))

save(myqtl1, myqtl2, file = "~/Dropbox/DO/R/myqtl1.rdt")

map5.dt <- data.frame(pos.gene = pos.gene[myqtl2.1[, 1]], pos.snp = pos.snp[myqtl2.1[, 2]], lod = myqtl2.1[, 3])

myMap(map5.dt)

mygene.sel <- rownames(myqtl1.1)
mygene.sel <- mygene.sel[ens.ucsc[mygene.sel, ]$chrom != 20]  # take off the X chromosome
mysnp.sel <- snp.ids[unique(myqtl2.1[, 2])]
mysnp.sel <- mysnp.sel[snp.inf[mysnp.sel, ]$Chr != 20]
myqtl.sel <- myqtl1.lod[which(gene.ids %in% mygene.sel), ]

par(mfrow = c(1, 3), lwd = 5, cex = 1.2)
venn(list(DOQTL = unique(qtl2.1[, 1]), MYQTL = unique(myqtl2.1[, 1])))
venn(list(DOQTL = unique(qtl2.1[, 2]), MYQTL = unique(myqtl2.1[, 2])))
venn(list(DOQTL = paste(qtl2.1[, 1], qtl2.1[, 2], sep = "-"), 
          MYQTL = paste(myqtl2.1[, 1], myqtl2.1[, 2], sep = "-")))

doqtl.sel.lod <- doqtl.lod[gene.sel1.ids, ]
myqtl.sel.lod <- myqtl1.lod[gene.sel1.ids, ]

hist.dt <- data.frame(lod = c(as.vector(doqtl.sel.lod), as.vector(myqtl.sel.lod)), 
                      group = c(rep("DOQTL", 288600), rep("MYQTL", 288600)))

ggplot(hist.dt, aes(x = lod, fill = group)) + geom_density(alpha = .3) + 
  theme_bw() + xlim(c(0, 15))
par(mfrow = c(1, 2))
for (i in c(2, 8)) {
  plot(doqtl.sel.lod[i, ], myqtl.sel.lod[i, ], pch = 20, cex = .2, xlab = "DOQTL", ylab = "MYQTL", main = i)
}
