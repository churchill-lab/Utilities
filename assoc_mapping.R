###########
# Libraries
library(qtl2)
library(qtl2db)
library(qtl2convert)
################################################################################
# Arguments:
# probs: genoprobs object for one chromosome in qtl2 format.
# pheno: matrix of phenotypes, samples in rows, phenotypes in columns.
# idx: column index of the phenotype that you want to map.
# addcovar: covariates matrix as used in scan1().
# intcovar: covariate to interact with QTL.
# K: kinship matrix for the current chromosome.
# markers: data.frame containing 4 columns: marker, chr, bp, cM.
# chr: Chromosome to map on.
# start: start position for mapping in Mb.
# end: end position for mapping in Mb.
# ncores: number of cores to use in mappin.
# snp.file: Location of the mySQL database containing the Sanger SNPs. This
#           was built by Karl and given to us.
assoc_mapping = function(probs, pheno, idx, addcovar, intcovar = NULL, K, 
                markers, chr, start, end, ncores = 1, 
                snp.file = "cc_variants.sqlite") {

  # Make sure that we have only one chromosome.
  if(length(probs) > 1) {
    stop(paste("Please provide a probs object with only the current chromosome."))
  } # if(length(probs) > 1)

  if(length(K) > 1) {
    stop(paste("Please provide a kinship object for the current chromosome."))
  } # if(length(K) > 1)

  stopifnot(length(probs) == length(K))

  # Create a function to query the SNPs.
  query_variants = create_variant_query_func(snp.file)

  # Convert marker positions to Mb if not already done.
  # The longest mouse chromosome is ~200 Mb, so 300 should cover it.
  if(max(markers[,3], na.rm = T) > 300) {
    markers[,3] = markers[,3] * 1e-6
  } # if(max(markers[,3], na.rm = T) > 300)

  # Split up markers into a vector of map positions.
  map = map_df_to_list(map = markers, pos_column = "pos")

  # Extract SNPs from the database
  snpinfo = query_variants(chr, start, end)

  # Names have to be replaced for future methods
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")

  # Index groups of similar SNPs.
  snpinfo = index_snps(map = map, snpinfo)

  # Keep samples that are not NA.
  keep = !is.na(pheno[,idx])

  # Convert genoprobs to snpprobs.
  snppr = genoprob_to_snpprob(probs[keep,], snpinfo)
  
  # Scan1.
  assoc = scan1(pheno = pheno[keep,idx, drop = FALSE], kinship = K[[1]][keep, keep],
          genoprobs = snppr, addcovar = addcovar[keep,,drop=FALSE], 
          intcovar = addcovar[keep,intcovar], cores = ncores)

  # Return the scan data.
  return(list(assoc, snpinfo))

} # assoc_mapping()

