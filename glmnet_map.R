################################################################################
# Perform association mapping in an interval, select SNPs above a certain LOD 
# and feed them into GLMnet to find groups of SNPs that may explain the QTL.
#
# Daniel Gatti
# dan.gatti@jax.org
# Aug. 4, 2017
################################################################################
library(qtl2)
library(glmnet)

source("assoc_mapping.R")

# Arguments:
# pheno: Numeric matrix containing the phenotype to map. Samples in rows and
#        1 column. Sample IDs must be in rownames.
# genoprobs: qtl2-style genoprobs object for all chromosomes.
# covar: Numeric matrix (as from model.matrix) containing the additive 
#        covariates to use in the model.
# map: List of markers positions for each chromosome.
# assoc: A scan1 object from assoc_mapping.
# snpinfo: A snpinfo object from assoc_mapping.
# lod.thr: Number indicating a LOD threshold to use when selecting SNPs.
#          We will only perform glmnet on SNPs with a LOD >= lod.thr.
# Returns: a list with 2 elements:
#          [[1]]: numeric matrix containing the shrunken coefficients;
#          [[2]]: snpinfo object that corresponds to the coeffients.
# The goal is to be able to use plot_snpasso() with the returned objects.
# NOTE: we are using alpha = 0.5, which is between LASSO and Ridge.
assoc_glmnet = function(pheno, genoprobs, covar, map, assoc, snpinfo, lod.thr = 1) {

  # Expand the SNPs and LOD scores.
  snpinfo.map = qtl2scan:::snpinfo_to_map(snpinfo)
  expand.lod  = qtl2plot:::expand_snp_results(assoc, snpinfo.map, snpinfo)
  snpinfo = data.frame(snpinfo, lod = expand.lod$lod, stringsAsFactors = FALSE)
  colnames(snpinfo)[ncol(snpinfo)] = "lod"

  # Select the SNPs with LOD scores above lod.thr.      
  print(paste(nrow(snpinfo), "SNPs input."))
  snpinfo = snpinfo[snpinfo$lod >= lod.thr,]
  print(paste(nrow(snpinfo), "SNPs with LOD > ", lod.thr, "."))

  # Get the SNP probs.
  new.snpinfo = index_snps(map = map, snpinfo = snpinfo)
  snpprobs = genoprob_to_snpprob(genoprobs, new.snpinfo)

  # Create an x matrix with the covariates and SNP probs.
  x = cbind(add_covar, snpprobs[[1]][,"A",])

  # Run glmnet and cross valdiation to select an optimal lambda.
  mod = glmnet(x = x, y = pheno[,1], alpha = 0.5)
  mod.cv = cv.glmnet(x = x, y = pheno[,1], alpha = 0.5)

  # Get the shrunken coefficients.
  coefs = as.matrix(coef(mod, s = mod.cv$lambda.min))

  # Strip off the additive coefficients.
  coefs = coefs[-(1:(ncol(covar)+1)),, drop = FALSE]

  return(list(coefs, new.snpinfo))

} # assoc_glmnet()


# Function to get 
get_nonzero_snps = function(glm_output) {

  coefs   = glm_output[[1]]
  snpinfo = glm_output[[2]]

  # Expand the SNPs and LOD scores.
  snpinfo.map = qtl2scan:::snpinfo_to_map(snpinfo)
  expand.snps = qtl2plot:::expand_snp_results(coefs, snpinfo.map, snpinfo)
  snpinfo = data.frame(snpinfo, expand.snps$lod, stringsAsFactors = FALSE)
  colnames(snpinfo)[ncol(snpinfo)] = "coef"

  # Select the SNPs with LOD scores above lod.thr.      
  print(paste(nrow(snpinfo), "SNPs input."))
  snpinfo = snpinfo[(snpinfo$coef) > 0,]
  print(paste(nrow(snpinfo), "SNPs with coefficients != 0."))

  return(snpinfo)

} # get_nonzero_snps()


#plot_snpasso(scan1output = abs(coefs), snpinfo = new.snpinfo, ylab = "coef")

#plot(mod, xvar = "dev", label = TRUE)
#plot(mod, xvar = "lambda", label = TRUE)

