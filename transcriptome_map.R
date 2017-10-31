################################################################################
# Given a set of eQTL, make a scatter plot of marker positions and gene 
# positions.
#
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 30, 2017
################################################################################
# Arguments: data: data.frame containing the following named columns:
#   chr_gene: character, the chromosome location of each transcript (or protein)
#   start_gene: numeric, the *megabase* start position of each transcript.
#   end_gene: numeric, the *megabase* end position of each transcript.
#   chr_qtl: character, the chromosome position of each QTL peak.
#   pos_qtl: numeric, the *megabase* position of each QTL peak.
#   lod: the LOD score of each peak.
#   cis: logical, TRUE if the peak is considered a local-eQTL.
tmap = function(data) {

  # Plot colors. First color is light and transparent, second is darker.
  trans.color = c(rgb(0,   0,   0, 0.3), rgb(  0,   0,  0,  0.9))
  cis.color   = c(rgb(0 ,0.3, 0.6, 0.3), rgb(0.1, 0.3, 0.7, 0.9))

  # Set up Chr lengths.
  data$chr_gene = factor(data$chr_gene , levels = c(1:19, "X", "Y", "MT"))
  data$chr_qtl  = factor(data$chr_qtl  , levels = c(1:19, "X", "Y", "MT"))
  chrlen = sapply(split(data$end_gene, data$chr_gene), max, na.rm = T)
  chrsum = c(0, cumsum(chrlen[-length(chrlen)]))
  names(chrsum) = names(chrlen)
  chrmid = chrsum + chrlen / 2

  # Create GMb positions.
  tgmb = (data$end_gene + data$start_gene) / 2
  qgmb = data$pos_qtl
  for(chr in levels(data$chr_gene)) {

    wh = which(data$chr_gene == chr)
    tgmb[wh] = tgmb[wh] + chrsum[chr]

    wh = which(data$chr_qtl == chr)
    qgmb[wh] = qgmb[wh] + chrsum[chr]

  } # for(i)

  # Get cis-qtl.
  col = rep(trans.color[1], nrow(data))
  col[data$cis == TRUE] = cis.color[1]

  data = cbind(data, qgmb, tgmb, col)
  data$col = as.character(data$col)

  trans = data[data$cis == FALSE,]
  cis   = data[data$cis == TRUE,]

  layout(matrix(1:2, 2, 1), height = c(0.7, 0.3))
  par(plt = c(0.09, 0.99, 0.05, 0.95))

  plot(trans$qgmb, trans$tgmb, col = 0, ann = F, axes = F, xaxs = "i", yaxs = "i")
  abline(h = chrsum, col = "grey80")
  abline(v = chrsum, col = "grey80")
  box()
  points(trans$qgmb, trans$tgmb, pch = 16, col = trans$col)
  points(cis$qgmb,   cis$tgmb,   pch = 16, col = cis$col)
  usr = par("usr")
  mtext(side = 1, line = 0.25, at = chrmid[1:20], text = names(chrmid)[1:20], 
        cex = 1.2)
  mtext(side = 2, line = 0.25, at = chrmid[1:20], text = names(chrmid)[1:20], 
        cex = 1.2, las = 2)
  mtext(side = 2, line = 2, text = "Transcript Position", cex = 1.5, las = 0)

  h.trans = hist(trans$qgmb, breaks = 0:ceiling(max(data$qgmb)), plot = F)
  h.cis   = hist(cis$qgmb,   breaks = 0:ceiling(max(data$qgmb)), plot = F)

  dens.trans = density(trans$qgmb, bw = 0.1, n = nrow(trans))
  dens.trans$y = nrow(trans) / sum(dens.trans$y) * dens.trans$y

  dens.cis = density(cis$qgmb, bw = 0.1, n = nrow(cis))
  dens.cis$y = nrow(cis) / sum(dens.cis$y) * dens.cis$y

  par(plt = c(0.09, 0.99, 0.15, 0.95))

  plot(dens.trans$x, dens.trans$y, type = "l", col = trans.color[2], xlim = usr[1:2],
       xaxt = "n", las = 2, xaxs = "i", ann = F)
  lines(dens.cis$x,  dens.cis$y,   col = cis.color[2])

  mtext(side = 1, line = 0.25, at = chrmid[1:20], text = names(chrmid)[1:20], 
        cex = 1.2)
  mtext(side = 1, line = 2, text = "QTL Position", cex = 1.5)
  mtext(side = 2, line = 2, text = "QTL Count", cex = 1.5, las = 0)

} # tmap()
