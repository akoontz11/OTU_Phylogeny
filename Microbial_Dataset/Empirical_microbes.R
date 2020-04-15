# Script for analysis of Bonnie's fungal dataset from CR

library(ape)
source("../MPDvsDelta.R")
setwd("/home/akoontz11/OTU_Phylogeny/Microbial_Dataset/")

# Read in and reformat fungal community matrix
fung.comm <- as.matrix(read.csv2("CRfungi_communitydat.csv", header=TRUE, sep=","))
soil.cores <- fung.comm[,1]
fung.comm <- apply(fung.comm[,-1],2,as.numeric)
rownames(fung.comm) <- soil.cores
colnames(fung.comm) <- gsub(pattern="X", replacement="", colnames(fung.comm))
sort(unique(c(fung.comm)))

# Read in fungal phylogenetic tree (aligned with mafft, built using RAxML)
fung.phylo <- read.tree("phylos/fungal/RAxML_bestTree.repCRfung_Feb28_phy")
# Prune tree of tips not present in the community matrix
fung.phylo <- keep.tip(fung.phylo, colnames(fung.comm))
plot(fung.phylo, show.tip.label=T)

# Find a ranking of sites by diversity before transforming the phylogeny
# (We use this as a vector to compare our mpd values against later)
# Create comparative comm object
fung.cc <- comparative.comm(fung.phylo, fung.comm, force.root = 0)
# Calculate mean pairwise distance value (for each of 24 sites in community)
fung.mpds <- .mpd(fung.cc, abundance.weighted=TRUE) 

# Build vector of delta values
deltas <- seq(0.1,3,by=0.1)
# Transform, building a matrix of mpds with sites as rows and delta values as columns
fung.test <- phy.d.transform.plot(fung.phylo, fung.comm, deltas)

# Summary of relative changes in diversity across sites
# Declare vector to receive average ranking changes
fung.siteShifts <- vector(length=length(deltas))
# Function for comparing site rankings at each delta value to original site rankings (prior to transformations)
compare.shifts <- function(x,d){
  # Generate a numeric vector of baseline site rankings (i.e. delta=1.0)
  baseline.rank <- rank(d)
  # Generate a vector to capture mean changes between rankings
  meanRankChanges <- vector(length = ncol(x))
  for (i in 1:ncol(x)){
    # Calculate the mean of the absolute changes between two different rankings, and place into vector
    meanRankChanges[i] <- mean(abs(baseline.rank - rank(x[,i])))
  }
  return(meanRankChanges)
}

fung.siteShifts <- compare.shifts(fung.test, fung.mpds)
fung.siteShifts
plot(fung.siteShifts ~ deltas, xlab="Delta values", ylab="Mean difference in site rankings from original", pch=16)
