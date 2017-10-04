################################################################################
# Simulate phenotypes with a given heritability.
# Ideas from Sweiger et al., Am J Hum Genet. 2016.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908190/
# 
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 3, 2017
################################################################################
library(mvtnorm)
library(qtl2)

# Retain only samples with data, i.e. no NA, Inf or Nan.
# All objects must have sample IDs in rownames.
# Arguments:
# pheno: numeric matrix containing a single column. Sample IDs must be in 
#        rownames.
# covar: numeric matrix of covariates as obtained from model.matrix(). Should
#        not include intercept. Sample IDs must be in rownames.
# K: numeric matrix containing kinship values as estimated by
#    qtl2geno::calc_kinship() or your favorite method. Samples IDs must be in 
#    rownames and colnames.
subset_data = function(pheno, covar, K) {

  samples = intersect(rownames(pheno), rownames(covar))
  samples = sort(intersect(rownames(K), samples))

  if(length(samples) == 0) {
    stop("No sample IDs in common across objects.")
  } # if(length(samples) == 0)

  pheno = pheno[samples,,drop = F]
  covar = covar[samples,,drop = F]
  K = K[samples, samples]

  tmp = cbind(pheno, covar, K)
  keep = complete.cases(tmp)
 
  if(length(keep) == 0) {
    stop("No non-missing data.")
  } # if(length(keep) == 0)

  pheno = pheno[keep,,drop = F]
  covar = covar[keep,,drop = F]
  K = K[keep, keep]

  return(list(pheno = pheno, covar = covar, K = K))

} # subset_data()

# Given phenotypes, covariates and a kinship matrix, estimate the 
# heritability and a confidence interval.
# All data must be complete and all rows should be aligned before calling
# this function.
# Arguments:
# pheno: numeric matrix containing a single column. Sample IDs must be in 
#        rownames.
# covar: numeric matrix of covariates as obtained from model.matrix(). Should
#        not include intercept. Sample IDs must be in rownames.
# K: numeric matrix containing kinship values as estimated by
#    qtl2geno::calc_kinship() or your favorite method. Samples IDs must be in 
#    rownames and colnames.
# nperm: integer denoting the number of simulations to perform.
# cores: integer denoting the number of cores to use for calculations.
h2_conf_int = function(pheno, covar, K, nperm = 1000, cores = 1) {

  stopifnot(rownames(pheno) == rownames(K))
  stopifnot(rownames(pheno) == rownames(covar))

  # Estimate the observed heritability.
  h2.obs = est_herit(pheno = pheno, kinship = K, addcovar = covar, reml = TRUE, 
           cores = cores)

  # Get the error covariance matrix by proportionally combining the kinship 
  # and identity matrices, based on the estimated heritability.
  n = nrow(pheno)
  twoK = 2 * K
  errcov = h2.obs * twoK + (1 - h2.obs) * diag(n)

  # Simulate nperm phenotypes and estimate their heritability.
  pheno.sim = t(rmvnorm(n = nperm, mean = rep(0, n), sigma = errcov))
  rownames(pheno.sim) = rownames(pheno)
  h2.perms = apply(pheno.sim, 2, est_herit, kinship = K, addcovar = covar, 
                         reml = TRUE, cores = cores)

  return(list(h2.obs = h2.obs, h2.perms = h2.perms))

} # h2_conf_int()


