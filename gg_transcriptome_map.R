################################################################################
# Transcriptome map using ggplot framework.
# Given a set of eQTL locations and (optionally) gene locations, plot the QTL
# location versus the gene location. Genes must be identified by Ensembl ID.
# The function can look up gene locations by Ensembl ID or the user may provide
# them.
# Points may be optionally colored by LOD score.
#
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 14, 2017
################################################################################
library(AnnotationHub)
library(tidyverse)

get_ensembl_genes = function() {

    hub = AnnotationHub()
    hub = query(hub, c("ensembl", "gtf", "mus musculus", "GRCm38"))
    hub = hub[grep("^Mus_musculus\\.GRCm38\\.[0-9]+\\.gtf$", hub$title)]
    latest.gtf = sort(hub$title)[length(hub)]
    ensembl = hub[[names(hub)[hub$title == latest.gtf]]]

    return(ensembl[ensembl$type == "gene"])

} # get_ensembl_genes()

# I had hoped to get the chr lengths from seqlengths(ensembl), but they
# don't always fill it in. I'm getting it from the maximum gene location
# on each chromosome.
get_chr_length = function(ensembl) {

  tmp = data.frame(chr = seqnames(ensembl), end = end(ensembl))
  tmp = tmp %>% filter(!(substring(chr, 1, 1) %in% c("G", "J"))) %>%
          group_by(chr) %>% 
          summarize(len = max(end))
  chrlen = tmp$len * 1e-6
  names(chrlen) = tmp$chr

  return(chrlen)

} # get_chr_length()


# Arguments:
# data: data.frame (or tibble) with the following columns:
#       ensembl: (required) character string containing the Ensembl gene ID.
#       qtl_chr: (required) character string containing QTL chromsome.
#       qtl_pos: (required) floating point number containing the QTL position 
#                in Mb.
#       qtl_lod: (optional) floating point number containing the LOD score.
#       gene_chr:  (optional) character string containing transcript chromsome.
#       gene_start: (optional) character string containing transcript start 
#                 postion in Mb.
#       gene_end:  (optional) character string containing transcript end
#                position in Mb.
# Returns:
# a plot of the QTL and gene location for each gene.
ggtmap = function(data, color.points = FALSE) {

  required.colnames = c("ensembl", "qtl_chr", "qtl_pos")
  if(all(!required.colnames %in% colnames(data))) {
    stop(paste("colnames must contain the following columns:", 
         paste(required.colnames, collapse = ",")))
  } # if(!required.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$ensembl = as.character(data$ensembl)
  data$qtl_chr = as.character(data$qtl_chr)

  # Get the latest Ensembl GTF.
  ensembl = get_ensembl_genes()
#  chr.len = get_chr_length(ensembl)

  gene.position.colnames = c("gene_chr", "gene_start", "gene_end")
  if(!all(gene.position.colnames %in% colnames(data))) {

    message(paste("Using Ensembl gene locations because optional gene",
            "position columns (", paste0(gene.position.colnames, collapse = ","),
            ") not found."))

    id    = ensembl$gene_id
    chr   = seqnames(ensembl)
    start = start(ensembl) * 1e-6
    end   = end(ensembl)   * 1e-6

    df = data.frame(ensembl = id, gene_chr = chr, gene_start = start, gene_end = end,
         stringsAsFactors = F)
    data = left_join(data, df, by = "ensembl")

  } # if(gene.position.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$gene_chr = as.character(data$gene_chr)

  # Get the gene mid-point.
  data = data %>% mutate(gene_pos = (gene_end + gene_start) * 0.5)

  # Fix the factor levels for the chr.
  all.chr = data %>% select(qtl_chr, gene_chr) %>%
            gather(k, v) %>%
            select(v) %>%
            distinct() %>%
            arrange(v)
  all.chr = all.chr$v[!is.na(all.chr$v)]

  if(length(grep("M", all.chr)) > 0) {
    wh = grep("M", all.chr)
    all.chr = all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
  }

  # Remove any NAs.
  data = data %>% na.omit

  data$qtl_chr  = factor(data$qtl_chr,  levels = all.chr[order(as.numeric(all.chr))])
  data$gene_chr = factor(data$gene_chr, levels = rev(all.chr[order(as.numeric(all.chr))]))


  print(ggplot(data, aes(x = qtl_pos, y = gene_pos)) +
    geom_point(aes(color = lod)) + {
      if(color.points) scale_color_continuous(low = "grey50", high = "red") 
    } +
    facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE) +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.text.x = element_text(angle = 90, hjust = 1)))

} # ggtmap()

# Test code
#data = read.csv("/hpcdata/gac/derived/Attie_DO_Islet_RNASeq/eQTL/additive/attie_DO378_additive_qtl_summary_thresh_6.csv")
#colnames(data) = c("ensembl", "marker", "qtl_chr", "qtl_pos", "lod")
#png("tmap.png", width = 2000, height = 2000, res = 164)
#ggtmap(data = data, color.points = TRUE)
#dev.off()
