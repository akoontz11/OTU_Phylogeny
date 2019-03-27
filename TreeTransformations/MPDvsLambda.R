# %%% PLOTTING PHYLOGENETIC DIVERSITY VALUES OF A GIVEN TREE OVER INCREASING VALUES OF PAGEL'S LAMBDA %%%

library(geiger)
library(pez)

# %%% USING LAJA DATASET %%%
data(laja)
# Creating a commmunity comparative ecology object based off of a phylogeny and a community matrix (species as columns, rows as communities)
data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
plot(invert.tree)

# "simulating evolution of discrete or continuous characters on a phylogenetic tree". 
# Returns "an array of simulated data, either two or three-dimensional. First dimension is the number of taxa, the second the number of characters, and the third the number of simulated data sets.
# Then, we're transposing the 2nd dimension
# Then, we're multiplying by 100 and rounding
# But...how does evolution of character traits simulate species abundances? Or perhaps this is just a means of obtaining those abundance values?
species.abundances <- round(t(sim.char(invert.tree, 5, model="BM", nsim=10)[,1,]) * 100)
# Remove negative abundances
species.abundances[species.abundances < 0] <- 0
# Make formatting pretty
colnames(species.abundances) <- invert.tree$tip.label
rownames(species.abundances) <- paste("site_", seq_len(nrow(species.abundances)))

# Function for capturing mpd values over series of lambda transformations
phy.l.transform <- function(phylo,l,nsim){
  # Create matrix with which to capture mpd.values vector for every loop iteration
  mpd.mat <- matrix(NA, nrow=nsim, ncol=length(l))
  colnames(mpd.mat) <- l
  rownames(mpd.mat) <- paste("Site",1:nsim,sep="_")
  for(i in 1:length(l)){
    # Apply phylogenetic transformation to tree
    s.phylo <- rescale(phylo, "lambda", l[i])
    # Create comparative data object
    c.data <- comparative.comm(s.phylo, species.abundances)
    # Calculate mean pairwise distance
    mpd.values <- .mpd(c.data, abundance.weighted=TRUE)
    # Store calculated mpd values within a matrix
    mpd.mat[,i] <- mpd.values
  }
  # For each 'site' or simulation, plot the change in mpd versus the increasing values of lambda
  # Specifying the range of mpd values
  ymin <- min(mpd.mat); ymax <- max(mpd.mat)
  # Plotting the first matrix row (i.e. Site1 values) 
  plot((mpd.mat[1,1:length(l)]) ~ l, xlab="lambda", ylab="MPD Values", ylim=c(ymin,ymax), pch=20)
  lines(l, mpd.mat[1,])
  # Below loop iterates through length of the matrix, adding connected points onto the plot.
  for(i in 2:nsim){
    points((mpd.mat[i,1:length(l)]) ~ l,col=i, pch=20)
    # Capture mpd values for current row of matrix, and connect points
    y <- (mpd.mat[i,])
    lines(l,y,col=i)
  }
  return(mpd.mat)
}
# Generating a vector of lambda values from 0 to 1
lambdas <- seq(0,1,by=0.05)
laja.test <- phy.l.transform(invert.tree, lambdas,10)

# Utilizing dataset from Barro Colorado Island (BCI) Tree Counts (in vegan package)
library(vegan)
data(BCI)
head(BCI)
class(BCI)
# For usage by pez, we need BCI to be a matrix
BCI <- as.matrix(BCI)

# BCI has rows as hectare plots, and columns as species names--this is our community matrix. 
# So, in order to create a comparative.comm object, we need a phylogenetic tree...

# Found this phylogeny online through the paper Kress et al. 2018, and while it comes from the same island, it doesn't have sites in common with the BCI dataset in vegan
tree_file <- "/home/austin/Documents/OTU_Phylogeny/TreeTransformations/BCI_Dataset/PhylogeneticResources/PhylogeneticResources/Vascular_Plants_rooted.dated.tre"
head(tree_file)
BCI.Tree <- read.tree(file=tree_file)
plot(BCI.Tree)

data <- comparative.comm(BCI.Tree, BCI)

# I'm doing something wrong here...it seems like this should be easier. Maybe there's a BCI phylogeny I can pull from somewhere in vegan? 
# Stuck...maybe ask Will how to get the required BCI data.

