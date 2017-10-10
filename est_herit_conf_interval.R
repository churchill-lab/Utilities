################################################################################
# Simulate phenotypes with a given heritability.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 3, 2017
################################################################################
library(mvtnorm)
library(qtl2)
library(ggplot2)

# Retain only samples with data, i.e. no NA, Inf or Nan.
# All objects must have sample IDs in rownames.
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

  h2.perms = est_herit(pheno.sim, kinship = K, addcovar = covar, 
                       reml = TRUE, cores = cores)

  return(list(obs = h2.obs, perms = h2.perms))

} # h2_conf_int()

# Function to plot the distribution and 95% confidence interval.
herit_dist_plot = function(h2, title = "") {

  mean.h2 = mean(h2$perms)
  q95 = quantile(h2$perms, probs = c(0.025, 0.975))

  # All data.
  df = data.frame(x = h2$perms)

  # Just the 95% confidence interval.
  dens = density(h2$perms)
  dens = data.frame(x = dens$x, y = dens$y)
  dens = subset(dens, x >= q95[1] & x <= q95[2])

  ggplot(data = df) + 
    geom_density(aes(x = x, y = ..density..), color = "black") +
    geom_area(data = dens, aes(x = x, y = y), 
    fill = rgb(1, 0, 0), alpha = 0.3) +
    geom_vline(xintercept = mean.h2) +
    labs(title = paste(title, "h2 =", format(mean.h2, digits = 3), "[",
         format(q95[1], digits = 3), ",", format(q95[2], digits = 3), "]"))

} # herit_dist_plot()


# From https://github.com/simecek/HPQTL/blob/master/R/normalize.matrix.r
# This uses the sample-variance method.
normalize.matrix = function(W) {

  stopifnot(class(W) == "matrix")
  stopifnot(is.numeric(W) & is.finite(W))

  n = nrow(W)
  k = 1 / n * (sum(diag(W)) - 2 * sum(W[upper.tri(W)])  / (n - 1))
  return(W / k)

} # normalize.matrix()

#load("/hpcdata/gac/projects/Svenson_DO850/data/svenson_DO850_qtl2_input_v2.Rdata")
#num.na = colSums(is.na(pheno.rz))
# Estimate K.
#K = calc_kinship(probs = genoprobs, type = "overall", cores = 4)
#rm(genoprobs)
#pheno = pheno[!is.na(pheno[,1]),, drop = F]
#covar = covar[rownames(pheno),]
#K = K[rownames(pheno), rownames(pheno)]


#tmp = subset_data(pheno = pheno.rz[,"bw.20",drop = F], covar = covar, K = K)
#h2.bw.20 = h2_conf_int(pheno = tmp$pheno, covar = tmp$covar, K = tmp$K, 
#                       nperm = 100, cores = 4)
