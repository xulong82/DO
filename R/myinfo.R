# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: show the results
# Rev: April 6, 2014
#---------------------------------------------------------------------------------------------------
# load("~/Dropbox/DO/R/myanalysis1.rdt")

chr.ucsc <- read.delim("~/Dropbox/Genome/chromInfo.txt", header = F)  # chromosomes length
chrlen <- cumsum(as.numeric(chr.ucsc$V2)) * 1e-6
names(chrlen) <- 1:21
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

ens.ucsc <- read.delim("~/Dropbox/Genome/myEnsGene.tx")
ens.ucsc <- ens.ucsc[!duplicated(ens.ucsc$name2), ]
ens.ucsc <- data.frame(row.names  = ens.ucsc$name2, 
                       chrom      = ens.ucsc$chrom,
                       chromStart = ens.ucsc$txStart,
                       chromEnd   = ens.ucsc$txEnd)
ens.ucsc$chrom <- as.character(ens.ucsc$chrom)
ens.ucsc <- ens.ucsc[-grep("random", ens.ucsc$chrom), ]
ens.ucsc <- ens.ucsc[-grep("chrUn", ens.ucsc$chrom), ]
ens.ucsc <- ens.ucsc[-grep("chrM", ens.ucsc$chrom), ]
ens.ucsc$chrom <- gsub("chrX", "chr20", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chrY", "chr21", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chr", "", ens.ucsc$chrom)
ens.ucsc$chrom <- as.numeric(ens.ucsc$chrom)

snp.ids <- data.frame(SNP_ID = colnames(founder.probs[1, , ]))
snp.inf <- merge(snp.ids, muga_snps, by = "SNP_ID", all.x = T, sort = F)
snp.inf$SNP_ID <- as.character(snp.inf$SNP_ID)
snp.inf$Chr <- as.character(snp.inf$Chr)
snp.inf$Chr <- gsub("X", "20", snp.inf$Chr)
snp.inf$Chr <- gsub("Y", "21", snp.inf$Chr)
snp.inf$Chr <- as.numeric(snp.inf$Chr)

chrlen1 <- c(0, chrlen)
ens.ucsc.pos <- chrlen1[ens.ucsc$chrom] + rowMeans(ens.ucsc[, 2:3]) * 1e-6
ens.ucsc <- cbind(ens.ucsc, pos = ens.ucsc.pos)

snp.inf.pos <- chrlen1[snp.inf$Chr] + snp.inf$Mb_NCBI38
snp.inf <- data.frame(row.names = snp.inf$SNP_ID, snp.inf[, -1], pos = snp.inf.pos)

save.image("~/Dropbox/DO/R/myinfo.rdt")