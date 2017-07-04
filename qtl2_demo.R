################################################################################
# qtl2 mapping example.
# Association mapping code from Petr Simecek.
# Daniel Gatti
# dan.gatti@jax.org
# June 2, 2017
################################################################################
options(stringsAsFactors = F)
library(qtl2)   # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
library(RSQLite)
library(dplyr)

#####################
# Load in the data. #
#####################
norm.dir = "/hpcdata/gac/projects/CGD_DO_Liver_RNASeq/normalized/"
expr = readRDS(paste0(norm.dir, "Svenson_DO474_normalized_genes_v3.rds"))
expr = t(expr)
annot = readRDS(paste0(norm.dir, "Svenson_DO474_gene_annotation_v3.rds"))

covar = readRDS(paste0(norm.dir, "Svenson_DO474_covar_batch_v3.rds"))
probs = readRDS("/hpcdata/gac/projects/Svenson_DO850/data/svenson_DO850_genoprobs_v2.rds")
probs = probs[rownames(expr),]
stopifnot(rownames(expr) == rownames(probs[[1]]))
markers = read.delim("/hpcdata/gac/derived/CGD_DO_Genoprobs/marker_grid_64K.txt")

# Not used here, but may prove handy.
# Convert DOQTL-style allele probabilities to qtl2-style data structure.
# genoprobs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = markers, 
#              chr_column = "chr", pos_column = "bp")
#
# Convert qtl2 36-state genotype probabilities to 8-state allele probabilities.
# qtl2geno::genoprob_to_alleleprob(probs = probs, cores = 8)

###################
# Kinship matrix. #
###################
K = calc_kinship(probs = probs, type = "loco", cores = 5)

######################
# Set up covariates. #
######################
addcovar = model.matrix(~Sex + Diet + Gen, data = covar)[,-1]

#########
# Scan1 #
#########
qtl = scan1(genoprobs = probs, pheno = expr[,1:10,drop = FALSE],
            kinship = K, addcovar = addcovar, cores = 10)

str(qtl)

################
# Plotting QTL #
################

# Create map.
map = split(markers$bp * 1e-6, markers$chr)
map = map[order(as.numeric(names(map)))]

# Linkage mapping plot.
plot_scan1(x = qtl, map = map, lodcolumn = 1, main = colnames(qtl)[1])

##################################################
# Calculate fouder allele effects (coefficients) #
##################################################

# Calculate linear model coefficients on Chr 3.
chr = 3
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno = expr[,1,drop = FALSE],
           kinship = K[[chr]], addcovar = addcovar)

# Plot only coefficients.
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     main = colnames(qtl)[1])

# Plot coefficients with LOD cruve.
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])

# Calculate BLUP coefficients on Chr 3.
qtl.blup = scan1blup(genoprobs = probs[,chr], pheno = expr[,1,drop = FALSE],
           kinship = K[[chr]], addcovar = addcovar, cores = 10)

plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])

########################
# Association mapping. #
########################
# NOTE: This still a bit messy. Improvements are coming...
# NOTE2: You need /data/gac/resource/CCsnps/ccfoundersnps.sqlite, which is 
#        a database of Sanger SNPs for the 8 DO founders.

# Arguments:
# probs: genoprobs object in qtl2 format.
# pheno: matrix of phenotypes, samples in rows, phenotypes in columns.
# idx: column index of the phenotype that you want to map.
# addcovar: covariates matrix as used in scan1().
# intcovar: covariate to interact with QTL.
# k: list of kinship matrices.
# markers: data.frame containing 4 columns: marker, chr, bp, cM.
# chr: Chromosome to map on.
# start: start position for mapping in Mb.
# end: end position for mapping in Mb.
# ncores: number of cores to use in mappin.
# db.file: Location of the mySQL database containing the Sanger SNPs.
assoc_mapping = function(probs, pheno, idx, addcovar, intcovar = NULL, k, markers, 
                chr, start, end, ncores, 
                db.file = "/data/gac/resource/CCsnps/ccfoundersnps.sqlite") {

  # Subset probs and K to keep only the current chromosome.
  probs = probs[,chr]
  k     = k[[chr]]

  # Split up markers into a vector of map positions.
  map = split(markers[,3] * 1e-6, markers[,2])
  nm  = split(markers[,1], markers[,2])
  map = mapply(function(x, y) { names(x) = y;x }, map, nm)
  map = map[order(as.numeric(names(map)))]

  # Extract SNPs from the database
  my_db = src_sqlite(db.file, create = FALSE)
  snpinfo = tbl(my_db, sql(paste0("SELECT * FROM snps WHERE chr='", 
            chr, "' AND pos_Mbp>=", start, " AND pos_Mbp<=", end))) %>%
            collect(n = Inf)

  # Names have to be replaced for future methods
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")

  # Index groups of similar SNPs.
  snpinfo = index_snps(map = map, snpinfo)

  # Find which phenotype data actually exist
  sel = !is.na(pheno[,idx])

  # Convert genoprobs to snpprobs.
  snppr = genoprob_to_snpprob(probs[sel,], snpinfo)
  
  # Scan1.
  assoc = scan1(pheno = pheno[sel,idx, drop = FALSE], kinship = k[sel,sel],
          genoprobs = snppr, addcovar = addcovar[sel,], 
          intcovar = addcovar[sel,intcovar], cores = ncores)

  # Return the scan data.
  return(list(assoc, snpinfo))

} # assoc_mapping()


assoc.3 = assoc_mapping(probs = probs, pheno = expr, idx = 1, 
          addcovar = addcovar, k = K, markers = markers, 
          chr = chr, start = 105, end = 110, ncores = 4)

# Association mapping plot.
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
            drop.hilit = 1)

# Association mapping plot with narrower boundaries.
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
            drop.hilit = 1, xlim = c(107, 109))

##############################################
# Adding a gene plot to association mapping. #
##############################################

# Create data structure for gene plot.
genes = data.frame(chr   = annot$chr,
                   start = annot$Start.Mbp * 1e6,
                   stop  = annot$End.Mbp * 1e6,
                   strand = annot$Strand,
                   Name   = annot$Symbol, stringsAsFactors = F)


layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
            drop.hilit = 1, xlim = c(107, 109))
# NOTE: If you don't subset the gene list, this takes a long time.
plot_genes(genes[genes$chr == chr & genes$start > 107e6 & genes$stop < 109e6,], 
           xlim = c(107, 109))

# Take control of the colors and plot the gene we mapped in red.
genes.ss = genes[genes$chr == chr & genes$start > 107e6 & genes$stop < 109e6,]
# NOTE: You have to sort the genes by position for this to work.
genes.ss = genes.ss[order(genes.ss$start),]
gene.colors = rep("black", nrow(genes.ss))
gene.colors[genes.ss$Name == "Gnai3"] = "red"

plot_genes(genes.ss, xlim = c(107, 109), colors = gene.colors)

plot_genes(genes[genes$chr == chr,], xlim = c(107, 109))


