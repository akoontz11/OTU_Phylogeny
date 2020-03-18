# %%% PLOTTING PHYLOGENETIC DIVERSITY VALUES OF A GIVEN TREE OVER INCREASING VALUES OF PAGEL'S DELTA %%%

library(geiger)
library(pez)

# Function for capturing mpd values over series of delta transformations
phy.d.transform <- function(phylo,abundance.matrix,d){
  # phylo--the phylogenetic tree to be transformed
  # d--the vector of delta values to utilize for the branch length transformation
  # abundance.matrix--species-sites abundance matrix, used to calculate mpd values
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # If `phylo` argument is NULL, return an mpd.mat object of NA
  if(is.null(phylo)){
    mpd.mat <- NA
    return(mpd.mat)
  } else{
    # Calculate number of sites from species-site abundance matrix
    nsim <- nrow(abundance.matrix)
    # Create matrix with which to capture mpd.values vector for every loop iteration
    mpd.mat <- matrix(NA, nrow=nsim, ncol=length(d))
    colnames(mpd.mat) <- d
    rownames(mpd.mat) <- paste("Site",1:nsim,sep="_")
    for(i in 1:length(d)){
      # Apply phylogenetic transformation to tree
      s.phylo <- rescale(phylo, "delta", d[i])
      # Create comparative data object
      # (Including the force.root argument, in order to handle unrooted phylogenies)
      c.data <- comparative.comm(s.phylo, abundance.matrix, force.root = 0)
      # Calculate mean pairwise distance
      mpd.values <- .mpd(c.data, abundance.weighted=TRUE) 
      # Store calculated mpd values within a matrix
      mpd.mat[,i] <- mpd.values
    }
    return(mpd.mat)
  }
}

# %%% Demonstration using laja dataset----
#data(laja)
#data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#plot(invert.tree)
# Updating species-site abundance matrix
#species.abundances <- round(t(sim.char(invert.tree, 5, model="BM", nsim=10)[,1,]) * 100)
#species.abundances[species.abundances < 0] <- 0
#colnames(species.abundances) <- invert.tree$tip.label
#rownames(species.abundances) <- paste("site_", seq_len(nrow(species.abundances)))
# Generating a vector of delta values from 0.1 to 3
#deltas <- seq(0.1,3,by=0.1)
# Demonstrating function
#laja.test <- phy.d.transform(invert.tree, deltas,species.abundances)

# # %%% Old version of function (with plotting commands)----
# # Function for capturing mpd values over series of delta transformations
# phy.d.transform.plot <- function(phylo,abundance.matrix,d){
#   # phylo--the phylogenetic tree to be transformed
#   # d--the vector of delta values to utilize for the branch length transformation
#   # abundance.matrix--species-sites abundance matrix, used to calculate mpd values
#   # title--title of plot generated (plotting mpd values against increasing delta values)
#   # If `phylo` argument is NULL, return an mpd.mat object of NA
#   if(is.null(phylo)){
#     mpd.mat <- NA
#     return(mpd.mat)
#   } else{
#     # Calculate number of sites from species-site abundance matrix
#     nsim <- nrow(abundance.matrix)
#     # Create matrix with which to capture mpd.values vector for every loop iteration
#     mpd.mat <- matrix(NA, nrow=nsim, ncol=length(d))
#     colnames(mpd.mat) <- d
#     rownames(mpd.mat) <- paste("Site",1:nsim,sep="_")
#     for(i in 1:length(d)){
#       # Apply phylogenetic transformation to tree
#       s.phylo <- rescale(phylo, "delta", d[i])
#       # Plot transformed phylogeny
#       plot(s.phylo, show.tip.label=FALSE, main=c("Delta = ",d[i]))
#       # Create comparative data object
#       # Including the force.root argument, in order to handle unrooted phylogenies
#       c.data <- comparative.comm(s.phylo, abundance.matrix, force.root = 0)
#       # Calculate mean pairwise distance
#       mpd.values <- .mpd(c.data, abundance.weighted=TRUE) # This means that abundance data in simulation will effect calculation of mpd...Does this makes sense, given what we're trying to do?
#       # Store calculated mpd values within a matrix
#       mpd.mat[,i] <- mpd.values
#     }
#     # PLOTTING--For each 'site' or simulation, plot the change in mpd versus the increasing values of delta
#     # Specifying the range values, removing any NAs
#     ymin <- min(mpd.mat, na.rm=T); ymax <- max(mpd.mat,na.rm=T)
#     # Plotting the first matrix row (i.e. Site1 values)
#     plot((mpd.mat[1,1:length(d)]) ~ d, xlab="delta", ylab="MPD Values", ylim=c(ymin,ymax), pch=20)
#     lines(d, mpd.mat[1,])
#     # Below loop iterates through length of the matrix, adding connected points onto the plot
#     for(i in 2:nsim){
#     points((mpd.mat[i,1:length(d)]) ~ d,col=i, pch=20)
#     # Capture mpd values for current row of matrix, and connect points
#     y <- (mpd.mat[i,])
#     lines(d,y,col=i)
#     }
#     return(mpd.mat)
#   }
# }