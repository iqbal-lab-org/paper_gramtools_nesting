if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

require(ape)
tr_DBL_MSP <- read.tree("/home/brice/Desktop/main_PhD/analyses/nesting_paper/analysis/outputs/msps_dimorphism/trees/DBL_DBLMSP/RAxML_bestTree.DBL_DBLMSP")
tr_MSP <- read.tree("/home/brice/Desktop/main_PhD/analyses/nesting_paper/analysis/outputs/msps_dimorphism/trees/DBLMSP/RAxML_bestTree.DBLMSP")

require(ggtree)
ggtree(tr_DBL_MSP) + geom_tiplab(size=3)
ggtree(tr_MSP) + geom_tiplab(size=3)
