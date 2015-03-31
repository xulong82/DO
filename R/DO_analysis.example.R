#A test DO analsis

library(capeDO)
library(capeAddons)

#===============================================================
#specify the directories in which each class of data sits
#===============================================================
base.data.dir <- "XXX" #a directory containing other data directories. This is where the data.obj will be saved after it is generated.
genotype.data.dir <- "XXX" #the directory where the genotype data are stored.
snp.data.dir <- "XXX" #the directory where the marker information table is stored
pheno.data.dir <- "XXX" #the directory where the phenotype data are stored
results.dir <- "XXX" #a directory where results will be stored
#===============================================================

#===============================================================
#specify parameters for the analysis
#===============================================================
phenotypes <- c("Sex", "Diet", "Weight1", "X..Fat1", "LTM1")
pheno.labels <- c("Weight", "% Fat", "Lean Mass")
pheno.to.covar <- c("Sex", "Diet")
single.covar <- c("Sex", "Diet")
eig.which <- c(1,2)
scan.what <- "eigentraits"
ref.allele <- "B"
allele.labels <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
use.pairs.thresh = TRUE
pairscan.thresh = 4
non.allelic.covar <- c("Sex", "Diet")
pair.perms = 2
max.pair.cor = 0.5
pval.correction <- "fdr"
p.or.q = 0.05
transform.to.phenospace <- TRUE
r2.thresh = 0.5


#===============================================================
# build the data object. This only needs to be done once
#===============================================================
#get the SNP information for loading into the genotype function
setwd(snp.data.dir)
marker.info = read.table("marker.info.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# build the genotype array from the HMM output files
geno <- read.geno(path = genotype.data.dir, marker.info = marker.info)
saveRDS(geno, "geno.RData")

pheno <- read.pheno(path = pheno.data.dir, filename = "svenson.phenotypes.v5.numeric.csv", delim = ",", na.strings = NA)

cross <- combine.geno.pheno(geno, pheno)

setwd(base.data.dir)
saveRDS(cross, "crossDO.RData")
#===============================================================


#===============================================================
# read in the cross and start an analysis
#===============================================================
setwd(results.dir)
cross <- readRDS("crossDO.RData")

setwd(results.dir)
cross <- select.pheno(cross, pheno.which = phenotypes)
cross <- create.covar(cross, pheno.which = pheno.to.covar)

layout.mat <- get.layout.mat(dim(cross$pheno)[2])
pdf("Pheno.by.Ind.pdf", width = (dim(layout.mat)[2]*4), height = (dim(layout.mat)[1]*4))
layout(layout.mat)
apply(cross$pheno, 2, plot)
dev.off()

pdf("Phenotype.Correlation.by.Sex.pdf")
plot.pheno.cor.panel(cross, color.by = "Sex", group.labels = c("Female", "Male"), text.cex = 1)
dev.off()

pdf("Phenotype.Correlation.by.Diet.pdf")
plot.pheno.cor.panel(cross, color.by = "Diet", group.labels = c("Chow", "HF"), text.cex = 1)
dev.off()

plot.pheno.dist(cross, pdf.label = "Phenotype.Distributions.pdf")
plot.pheno.qq(cross, pdf.label = "Phenotype.QQ.plots.Normalized.pdf")

cross <- norm.pheno(cross, mean.center = TRUE)

pdf("Phenotype.Correlation.by.Sex.Normalized.pdf")
plot.pheno.cor.panel(cross, color.by = "Sex", group.labels = c("Female", "Male"), text.cex = 1.5, pheno.labels = pheno.labels)
dev.off()

pdf("Phenotype.Correlation.by.Diet.Normalized.pdf")
plot.pheno.cor.panel(cross, color.by = "Diet", group.labels = c("Chow", "HF"), text.cex = 1.5, pheno.labels = pheno.labels)
dev.off()

plot.pheno.dist(cross, pdf.label = "Phenotype.Distributions.Normalized.pdf")
plot.pheno.qq(cross, pdf.label = "Phenotype.QQ.plots.Normalized.pdf")

cross <- get.eigentraits(cross, scale.pheno = FALSE, normalize.pheno = FALSE)


pdf("svd.pdf")
plotSVD(cross, orientation = "vertical")
dev.off()

cross <- select.eigentraits(cross, traits.which = eig.which)

cross <- singlescan(cross, n.perm = 2, ref.allele = ref.allele, 
                    covar = single.covar, alpha.for.pairs = 0.01, alpha.for.covar = 0.05, 
                    scan.what = scan.what, auto.covar.selection = FALSE, verbose = TRUE)

saveRDS(cross, "crossDO.RData")

pdf("Singlescan.Standardized.pdf", width = 20, height = 8)
plotSinglescan(cross, standardized = TRUE, view = "overview", allele.labels = allele.labels)
dev.off()

pdf("Singlescan.pdf", width = 20, height = 8)
plotSinglescan(cross, standardized = FALSE, view = "overview", include.covars = TRUE, allele.labels = allele.labels)
dev.off()

pdf("Singlescan.Standardized.no.Covar.pdf", width = 20, height = 8)
plotSinglescan(cross, standardized = TRUE, view = "overview", allele.labels = allele.labels, include.covars = FALSE)
dev.off()

pdf("Singlescan.no.Covar.pdf", width = 20, height = 8)
plotSinglescan(cross, standardized = FALSE, view = "overview", include.covars = FALSE, allele.labels = allele.labels)
dev.off()

cross <- set.pairscan.thresh(cross, pairscan.thresh = pairscan.thresh)

cross <- select.markers.for.pairscan(cross, use.pairs.threshold = use.pairs.thresh, pairscan.thresh = pairscan.thresh)

cross <- pairscan(cross, scan.what = scan.what, n.perm = pair.perms, max.pair.cor = max.pair.cor, verbose = TRUE, num.pairs.limit = 300000)

saveRDS(cross, "crossDO.RData")

plotPairscan(cross, phenotype = NULL, pdf.label = "Pairscan.Regression.pdf", verbose = TRUE)

cross <- error.prop(cross, perm = FALSE, verbose = TRUE)
cross <- error.prop(cross, perm = TRUE, verbose = TRUE)

saveRDS(cross, "crossDO.RData")

cross <- calc.p(cross, pval.correction = pval.correction)

cross <- direct.influence(cross, transform.to.phenospace = transform.to.phenospace, verbose = TRUE, pval.correction = pval.correction, save.permutations = TRUE)

saveRDS(cross, "crossDO.RData")

pdf("variant.influences.pdf")
plotVariantInfluences(cross, p.or.q = p.or.q, standardize = FALSE, not.tested.col = "lightgray")
dev.off()


writeVariantInfluences(cross, p.or.q = max(c(p.or.q, 0.2)), filename = "Variant.Influences.csv")


cross <- get.network(cross, p.or.q = p.or.q, collapse.linked.markers = FALSE)
cross <- get.network(cross, p.or.q = p.or.q, r2.thresh = r2.thresh)

saveRDS(cross, "crossDO.RData")

pdf(paste("Network.Collapsed.r2.",r2.thresh, ".pdf", sep = ""), width = 12, height = 7)
plotNetwork(cross)
dev.off()


pdf(paste("VariantInfluences.Collapsed.Net.r2.", r2.thresh, ".from.pdf", sep = ""))
plotCollapsedVarInf(cross, expand.labels = FALSE, all.markers = FALSE)
dev.off()


node.size <- arrow.offset <- 1.5; label.offset <- node.size*1.2
legend.radius <- 9; label.cex = 1
pdf("Network.View.Full.pdf", width = 40, height = 40)
plotNetwork2(cross, p.or.q = p.or.q, collapsed.net = FALSE, r2.thresh = r2.thresh, node.radius = node.size, label.nodes = TRUE, label.offset = label.offset, label.cex = label.cex, legend.radius = legend.radius, legend.cex = 1.8, arrow.offset = arrow.offset)
dev.off()


node.size <- arrow.offset <- 4; edge.lwd = 3
label.offset <- node.size*1.2
legend.cex = 2; label.cex = 2
legend.radius <- 10; net.layout = NULL
pdf("Network.View.Collapsed.pdf", width = 20, height = 20)
net.results <- plotNetwork2(cross, p.or.q = p.or.q, collapsed.net = TRUE, r2.thresh = r2.thresh, node.radius = node.size, layout.matrix = net.layout, label.nodes = TRUE, label.offset = label.offset, label.cex = label.cex, legend.radius = legend.radius, legend.cex = legend.cex, arrow.offset = arrow.offset, edge.lwd = edge.lwd)
dev.off()

plotInteractions(cross)


#====================================================================================================================================
#====================================================================================================================================
