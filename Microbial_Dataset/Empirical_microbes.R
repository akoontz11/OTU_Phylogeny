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
fung.ses.mpds <- .ses.mpd(fung.cc, abundance.weighted=TRUE) 
fung.ses.mpds <- fung.ses.mpds$mpd.obs.z

# Build vector of delta values
deltas <- seq(0.1,3,by=0.1)
# Transform, building a matrix of mpds with sites as rows and delta values as columns
fung.test <- phy.d.transform.plot(fung.phylo, fung.comm, deltas, plot.title="")

# %%% SESmpd CORRELATIONS %%%----
# Function for measuring correlation between original SESmpd value and SESmpd values for each site 
SESmpd.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}
# Measure correlations for each delta value
fung.correl <- SESmpd.correls(fung.test, fung.ses.mpds)
fung.correl
plot(fung.correl ~ deltas, xlab="Delta values", ylab="Correlations to original SESmpd", pch=16)

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

fung.siteShifts <- compare.shifts(fung.test, fung.ses.mpds)
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
laja.ses.mpds <- .ses.mpd(laja.cc, abundance.weighted=TRUE) 
laja.ses.mpds <- laja.ses.mpds$mpd.obs.z

# Generating a vector of delta values from 0.1 to 3
deltas <- seq(0.1,3,by=0.1)
# Demonstrating function
laja.test <- phy.d.transform.plot(invert.tree, species.abundances, deltas, plot.title="")

# %%% MPD CORRELATIONS %%%----
# Function for measuring correlation between original SESmpd value and SESmpd values for each site 
SESmpd.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}
# Measure correlations for each delta value
laja.correl <- SESmpd.correls(laja.test, laja.ses.mpds)
laja.correl
plot(laja.correl ~ deltas, xlab="Delta values", ylab="Correlations to original SESmpd", pch=16)

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

laja.siteShifts <- compare.shifts(laja.test, laja.ses.mpds)
laja.siteShifts
plot(laja.siteShifts ~ deltas, xlab="Delta", ylab="Mean rank shifts", pch=16, )

# %%% PLOTTING MPD FOR BOTH DATASETS %%%----
library(RColorBrewer)
library(viridis)
# Code derived from plotting code for phy.d.transform.plot, in MPDvsDelta.R

# Raw values
# Plotting parameters: two rows, one column
par(mfcol = c(2, 1), mar = c(0,0,0,0), oma = c(4, 4, .5, .5), 
    mgp = c(2, 0.6, 0))
# CR fungi
ymin <- min(fung.test, na.rm=T)-0.5; ymax <- max(fung.test,na.rm=T)+0.5
fung.colors <- brewer.pal(9, "YlOrRd")
# plot((fung.test[1,1:length(deltas)]) ~ deltas, ylim=c(ymin,ymax), pch=20, axes=FALSE,
#      xlab="Delta", ylab="SESmpd", main="Costa Rican fungal communities")
plot((fung.test[1,1:length(deltas)]) ~ deltas, ylim=c(ymin,ymax), pch=20, axes=FALSE)
lines(deltas, fung.test[1,])
# Below loop iterates through length of the matrix, adding connected points onto the plot
for(i in 2:nrow(fung.test)){
  points((fung.test[i,1:length(deltas)]) ~ deltas,col=fung.colors[i], pch=20)
  # Capture mpd values for current row of matrix, and connect points
  y <- (fung.test[i,])
  lines(deltas,y,col=fung.colors[i])
}
axis(1L, labels = FALSE, tck=-0.02)
axis(1L, labels = FALSE, tck=0.02)
axis(2L)
box()
# Laja
ymin <- min(laja.test, na.rm=T)-0.5; ymax <- max(laja.test,na.rm=T)+0.5
laja.colors <- brewer.pal(9, "Greens")
# plot((laja.test[1,1:length(deltas)]) ~ deltas, ylim=c(ymin,ymax), pch=20, axes=FALSE, 
#      xlab="Delta", ylab="SESmpd", main="Laja")
plot((laja.test[1,1:length(deltas)]) ~ deltas, ylim=c(ymin,ymax), pch=20, axes=FALSE)
lines(deltas, laja.test[1,])
# Below loop iterates through length of the matrix, adding connected points onto the plot
for(i in 2:nrow(laja.test)){
  points((laja.test[i,]) ~ deltas, col=laja.colors[i], pch=20)
  # Capture mpd values for current row of matrix, and connect points
  y <- (laja.test[i,])
  lines(deltas,y,col=laja.colors[i])
}
axis(1L, tck=-0.02)
axis(1L, tck=0.02)
axis(2L)
box()
mtext("Delta", side = 1, outer = TRUE, line = 2.2, cex = 1.5)
mtext("SESmpd", side = 2, outer = TRUE, line = 2.2, cex = 1.5)

# Plotting rank shifts per site
# Plotting parameters: two rows, one column
par(mfcol = c(2, 1), mar = c(0,0,0,0), oma = c(4, 4, .5, 0.5), 
    mgp = c(2, 0.6, 0), xpd=FALSE)
# CR fungi
fung.siteShifts <- compare.shifts(fung.test, fung.ses.mpds)
ymin <- min(fung.siteShifts)-0.5; ymax <- max(fung.siteShifts)+0.5
plot(fung.siteShifts ~ deltas, axes=FALSE, pch=20, col=fung.colors[6], ylim=c(ymin,ymax))
axis(1L, labels = FALSE, tck=-0.02)
axis(1L, labels = FALSE, tck=0.02)
axis(2L)
box()
# text(x = 2.2, y = 2.0, "Costa Rica \n fungal data", cex=0.8)
mtext("Costa Rica \n fungal data", side = 1, outer = FALSE, line = -1.0, adj=0, cex = 0.95, at=c(6.5,2.1))
# Laja
laja.siteShifts <- compare.shifts(laja.test, laja.ses.mpds)
ymin <- min(laja.siteShifts)-0.5; ymax <- max(laja.siteShifts)+0.5
plot(laja.siteShifts ~ deltas, axes=FALSE, pch=20, col=laja.colors[6], ylim=c(ymin,ymax))
axis(1L, tck=-0.02)
axis(1L, tck=0.02)
axis(2L)
box()
# text(x = 2.2, y = 2.0, "Laja data", cex=0.8)
mtext("Laja \n macroinvertebrate data", side = 1, outer = FALSE, line = -1.0, adj=0, cex = 0.95, at=c(6.5,2.1))
mtext("Delta", side = 1, outer = TRUE, line = 2.2, cex = 1.3)
mtext("Mean rank shifts", side = 2, outer = TRUE, line = 2.2, cex = 1.3)
