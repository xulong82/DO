################################################################################
# Create a transcriptome map, given a set of eQTL with transcript locations.
# This treats the QTL and the transcript as points, not intervals.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 19, 2012
################################################################################
library(org.Mm.eg.db)

# This function currently works only for the mouse.
# We require a large matrix of QTL for each transcript with transcript
# identifiers of some kind.  If no transcript identifiers are provided, we
# attempt to get the transcript locations from package org.Mm.eg.db.  If this
# fails, an error is returned.  Optionally, if you have annotation for your
# transcripts, you may provide it in the annot argument and we will match
# trancripts found in the QTL matrix with those in the annotation matrix.
# Arguments: qtl: data.frame, with the following columns:
#                 1: Transcript ID: character, identifier for transcripts, could 
#                    be Entrez Gene, Ensemble, MGI, gene symbols, etc.
#                 2: SNP ID:, character, the ID of the SNP where the QTL occurs.
#                 3: chr: character, the Chr on which th eQTL lies.
#                 4: pos: numeric, the position of the QTL on the chr in bp.
#                 5: score: numeric, the significance statistic for the QTL. The
#                    names of the statisic will be taken from the name of this
#                    column (ie, LOD, p-value, etc.)
#                 If the annot argument is empty, you *must* include the 
#                 annot.type argument to indicate the type of the transcript IDs.
#            annot: (optional), data.frame with transcript annotation. Must have
#                   the following columns:
#                 1: Gene ID: character, identifier for transcripts, could be Entrez
#                    transcript, Ensemble, MGI, transcript symbols, etc. Must be the same as
#                    the IDs used for the QTL.
#                 2: chr: character, the chromosome on which the transcript lies
#                    or NA if not known.
#                 3: start: numeric, the start of the transcript in bp.
#                 4: end: numeric, the end of the transcript in bp.
#                 If you provide the annot argument, then you will not need the
#                 annot.type argument.
#            annot.type: (optional), character, one of (ensembl, entrez, mgi or 
#                        refseq).  The type of annotation provided in the qtl
#                        argument.
#                 If the annot argument is empty, you *must* include the 
#                 annot.type argument to indicate the type of the transcript IDs.
#            ...: additional arguments to pass onto plot.
transcriptome.map = function(qtl, annot, annot.type = c("ensembl", "entrez", 
                             "mgi", "refseq"), ...) {
  # QTL not NULL.
  if(is.null(qtl)) {
    stop(paste("qtl is empty.  You must provide some QTL data."))
  } # if(is.null(qtl))

  # QTL and annot.type or QTL and annot provides.
  if(is.null(annot.type) && is.null(annot)) {
    stop(paste("You must either provide the annotation data frame or the",
         "type of annotation in annot.type."))
  } # if(is.null(annot.type) & is.null(annot))

  # If we have annotation from the user, then use it.
  gmb = NULL
  if(!is.null(annot)) {

    # See if we have transcript locations for all of the QTL transcripts.
    missing = which(!qtl[,1] %in% annot[,1])
    if(length(missing > 0)) {
      warning(paste("There are", length(missing), "QTL with no transcript",
              "location information.  Removing these from the plot."))
      qtl = qtl[-missing,]
    } # if(length(missing > 0))

    # Change "X" to "20" if we see it.
    qtl[,3]   = sub("X", "20", qtl[,3])
    annot[,2] = sub("X", "20", annot[,2])

    # Match the QTL transcripts with their locations.
    m = match(qtl[,1], annot[,1])

    # Append the trancript locations onto the QTL matrix.
    qtl = cbind(qtl, annot[m,])

    # Use the mean transcript location.
    qtl = cbind(qtl, mean = rowMeans(qtl[,8:9]))

    # Create a trancript and marker location matrix.
    gmb = matrix(c(qtl[,4], qtl[,10]) * 1e-6, nrow = nrow(qtl), ncol = 2,
          dimnames = list(qtl[,1], c("Marker (GMb)", "Transcript (GMb)")))
    gmb = cbind(gmb, qtl[,5])
    colnames(gmb)[3] = colnames(qtl)[5]

  } else {
  # Otherwise, use org.Mm.eg to try to extract annotation.
    # DMG: To be added later...
  } # else

  # Get the mouse Chr lengths.
  chrlen = org.Mm.egCHRLENGTHS
  chrlen = chrlen[-grep("random", names(chrlen))]
  chrlen = cumsum(chrlen * 1e-6)
  chrlen = chrlen[names(chrlen) %in% c(1:19, "X")]
  names(chrlen) = sub("X", "20", names(chrlen))
  chrmid = diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

  # Place the genome Mb locations in the gmb matrix.
  for(c in 2:20) {
    rng = which(qtl[,3] == c)
    gmb[rng,1] = gmb[rng,1] + chrlen[c-1]
    rng = which(qtl[,7] == c)
    gmb[rng,2] = gmb[rng,2] + chrlen[c-1]
  } # for(c in 2:20)

  # Change the "20" in chrlen names back to "X".
  names(chrlen) = sub("20", "X", names(chrlen))
  
  # Create the plot.
  par(font = 2, font.lab = 2, font.axis = 2, las = 1,
      plt = c(0.12, 0.95, 0.12, 0.9), lend = 2)
  plot(gmb[,1:2], xlim = c(0, max(chrlen)), ylim = c(0, max(chrlen)), col = 0,
       xaxs = "i", yaxs = "i")
  abline(h = chrlen, col = "grey80")
  abline(v = chrlen, col = "grey80")
  usr = par("usr")
  par(xpd = NA)
  txt.loc = usr[2] + (usr[2] - usr[1]) * 0.02
  text(chrmid, txt.loc, names(chrlen))
  text(txt.loc, chrmid, names(chrlen))
  par(xpd = F)
  breaks = quantile(gmb[,3], 0:100/100)
  rcol = grey(0:100 / 100)
#  rcol = colorRampPalette(c("cyan", "blue"))(100)
#  rcol = heat.colors(100)
  plot.cols = rep(rcol[1], nrow(qtl))
  for(i in 1:100) {
    plot.cols[gmb[,3] > breaks[i]] = rcol[i]
  } # for(i)
  points(gmb[,1:2], pch = 20, col = plot.cols, cex = 0.7)
  par(xpd = NA)
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2, border = 1)
  
  # Create a legend.
#  par(xpd = NA)
#  legend(-400, -175, format(quantile(qtl[,4]), digits = 2), pch = 16, col = 
#         rcol[c(1, 25, 50, 75, 100)], horiz = T)
#  rect(-400, -400, 650, -175, col = "white", border = 1)
#  x.pos = c(-250, -50, 150, 350, 550)
#  points(x.pos, rep(-225, 5), pch = 16, col = rcol[c(1, 25, 50, 75, 100)])
#  text(x.pos, rep(-300, 5), format(quantile(qtl[,4]), digits = 2), cex = 0.8)
  par(xpd = F)

} # transcriptome.map()


################################################################################
# QTL histogram that counts the number of QTL within a siding window.
# Daniel Gatti
# Dan.Gatti@jax.org
# Mar. 21, 2012
# Arguments: qtl: data.frame, with the following columns:
#                 1: Transcript ID: character, identifier for transcripts, could 
#                    be Entrez Gene, Ensemble, MGI, gene symbols, etc.
#                 2: SNP ID:, character, the ID of the SNP where the QTL occurs.
#                 3: chr: character, the Chr on which th eQTL lies.
#                 4: pos: numeric, the position of the QTL on the chr in bp.
#                 5: score: numeric, the significance statistic for the QTL. The
#                    names of the statisic will be taken from the name of this
#                    column (ie, LOD, p-value, etc.)
#             width: numeric, window width in Mb. Default = 1.
#             step: numeric, slide length in Mb. Default = 1.
################################################################################
qtl.hist = function(qtl, width = 1, step = 1) {

  # Get the mouse Chr lengths.
  chrlen = org.Mm.egCHRLENGTHS
  chrlen = chrlen[-grep("random", names(chrlen))]
  chrlen = chrlen[names(chrlen) != "Y" & names(chrlen) != "M"]
  chrlen = chrlen * 1e-6  
  chrsum = cumsum(chrlen)
  chrsum = chrsum[names(chrsum) %in% c(1:19, "X")]
  names(chrsum) = sub("X", "20", names(chrsum))
  chrmid = diff(c(0, chrsum)) * 0.5 + c(0, chrsum[-length(chrsum)])

  # Convert positions to Mb, if required.
  if(max(qtl[,4]) > 200) {
    qtl[,4] = qtl[,4] * 1e-6
  } # if(max(qtl[,4] > 200)

  # Split the QTL into chromosomes.
  qtl = split(qtl, qtl[,3])
  qtl = qtl[order(as.numeric(names(qtl)))]

  h = NULL
  for(i in 1:length(qtl)) {
    start = seq(0, (chrlen[[i]] - width), step)
    end   = start + width
    mids = rowMeans(cbind(start, end))
    mat = outer(qtl[[i]][,4], start, ">=") & outer(qtl[[i]][,4], end, "<")
    h = rbind(h, cbind(rep(i, length(mids)), mids, colSums(mat)))
  } # for(i)

  # Add the chromosome lengths to convert to GMb.
  h[,2] = h[,2] + c(0, chrsum)[h[,1]]

  # Make the plot.
  par(font = 2, font.axis = 2, font.lab = 2, las = 1)
  plot(h[,2], h[,3], col = 0, xaxs = "i", xlab = "Genome Mb",
       ylab = "Distant QTLs")
  abline(v = chrsum, col = "grey")
  points(h[,2], h[,3], type = "h")
  usr = par("usr")
  text(chrmid, usr[4] * 0.92, names(chrlen), cex=0.8, col="grey50")

  colnames(h) = c("Chr", "Midpoint", "Num.eQTL")
  return(h)
} # qtl.hist()
