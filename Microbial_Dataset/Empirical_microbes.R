# Script for analysis of empirical datasets: Bonnie's fungal dataset from CR, and laja dataset

library(ape)
setwd("/home/akoontz11/OTU_Phylogeny/Microbial_Dataset/")
source("../MPDvsDelta.R")

# --- %%%% FUNGAL DATASET %%%%----
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
fung.test <- phy.d.transform.plot(fung.phylo, fung.comm, deltas, plot.title="Shift in diversity of fungal communities over Delta values")

# %%% MPD CORRELATIONS %%%----
# Function for measuring correlation between original mpd value and mpd values for each site 
mpd.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}
# Measure correlations for each delta value
fung.correl <- mpd.correls(fung.test, fung.mpds)
fung.correl
plot(fung.correl ~ deltas, xlab="Delta values", ylab="Correlations to original mpd", pch=16)

# %%% SHIFTS IN SITE RANKINGS %%%----
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

# ----%%%% LAJA DATASET %%%%----
data(laja)
plot(invert.tree)

# Updating species-site abundance matrix (rather than using river.sites variable)
species.abundances <- round(t(sim.char(invert.tree, 5, model="BM", nsim=10)[,1,]) * 100)
species.abundances[species.abundances < 0] <- 0
colnames(species.abundances) <- invert.tree$tip.label
rownames(species.abundances) <- paste("site_", seq_len(nrow(species.abundances)))

# Create comparative comm object
laja.cc <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
# Calculate mean pairwise distance value
laja.mpds <- .mpd(laja.cc, abundance.weighted=TRUE) 

# Generating a vector of delta values from 0.1 to 3
deltas <- seq(0.1,3,by=0.1)
# Demonstrating function
# laja.test <- phy.d.transform.plot(invert.tree, river.sites, deltas, plot.title="Shift in diversity of Rio Laja communities over delta values")
laja.test <- phy.d.transform.plot(invert.tree, species.abundances, deltas, plot.title="Shift in diversity of Rio Laja macroinvertebrate communities over Delta values")

# %%% MPD CORRELATIONS %%%----
# Function for measuring correlation between original mpd value and mpd values for each site 
mpd.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}
# Measure correlations for each delta value
laja.correl <- mpd.correls(laja.test, laja.mpds)
laja.correl
plot(laja.correl ~ deltas, xlab="Delta values", ylab="Correlations to original mpd", pch=16)

# %%% SHIFTS IN SITE RANKINGS %%%----
# Declare vector to receive average ranking changes
laja.siteShifts <- vector(length=length(deltas))
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

laja.siteShifts <- compare.shifts(laja.test, laja.mpds)
laja.siteShifts
plot(laja.siteShifts ~ deltas, xlab="Delta values", ylab="Mean difference in site rankings from original", pch=16)