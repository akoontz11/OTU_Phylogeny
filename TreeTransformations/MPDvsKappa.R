# %%% PLOTTING PHYLOGENETIC DIVERSITY VALUES OF A GIVEN TREE OVER INCREASING VALUES OF PAGEL'S KAPPA %%%

library(geiger)
library(pez)

# %%% USING LAJA DATASET %%%
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
plot(invert.tree)

species.abundances <- round(t(sim.char(invert.tree, 5, model="BM", nsim=10)[,1,]) * 100)
species.abundances[species.abundances < 0] <- 0
colnames(species.abundances) <- invert.tree$tip.label
rownames(species.abundances) <- paste("site_", seq_len(nrow(species.abundances)))

# Function for capturing mpd values over series of kappa transformations
phy.k.transform <- function(phylo,k,nsim){
  # Create matrix with which to capture mpd.values vector for every loop iteration
  mpd.mat <- matrix(NA, nrow=nsim, ncol=length(k))
  colnames(mpd.mat) <- k
  rownames(mpd.mat) <- paste("Site",1:nsim,sep="_")
  for(i in 1:length(k)){
    # Apply phylogenetic transformation to tree
    s.phylo <- rescale(phylo, "kappa", k[i])
    # Create comparative data object
    c.data <- comparative.comm(s.phylo, species.abundances)
    # Calculate mean pairwise distance
    mpd.values <- .mpd(c.data, abundance.weighted=TRUE)
    # Store calculated mpd values within a matrix
    mpd.mat[,i] <- mpd.values
  }
  # For each 'site' or simulation, plot the change in mpd versus the increasing values of kappa
  # Specifying the range values
  ymin <- min(mpd.mat); ymax <- max(mpd.mat)
  # Plotting the first matrix row (i.e. Site1 values) 
  plot((mpd.mat[1,1:length(k)]) ~ k, xlab="kappa", ylab="MPD Values", ylim=c(ymin,ymax), pch=20)
  lines(k, mpd.mat[1,])
  # Below loop iterates through length of the matrix, adding connected points onto the plot.
  for(i in 2:nsim){
    points((mpd.mat[i,1:length(k)]) ~ k,col=i, pch=20)
    # Capture mpd values for current row of matrix, and connect points
    y <- (mpd.mat[i,])
    lines(k,y,col=i)
  }
  return(mpd.mat)
}
# Generating a vector of kappa values from 0 to 3 with 0.1 increments
kappas <- seq(0,3,by=0.1)
kappas

laja.test <- phy.k.transform(invert.tree, kappas,10)
