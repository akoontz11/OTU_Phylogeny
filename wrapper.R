library(pez)
library(ape)
library(geiger)
library(phytools)
library(parallel)

# ---------------------------
# %%% ADD BRANCH FUNCTION %%%
# ---------------------------
# Function for appending trees onto the tips of an existing tree. 
# The birth, death, and steps arguments all describe the parameters of the simulated trees to be appended to the original tree.
add.branch <- function(tree,birth,death,steps,type=c("pops","seq.err")){
  # Set type of branch addition according to argument (this specifies only the letters used for the tips of the appended tree)
  type <- match.arg(type)
  # Determine number of species present in original tree. Additional trees will be added to each originally present tip.
  tips <- Ntip(tree)
  # Make copy of original tree. This is in order to keep track of the tip names on the original tree.
  orig.tree <- tree
  # For each species in original tree,
  for(i in 1:tips){
    # get name of original tip.
    orig.tip <- orig.tree$tip.label[i]
    # Generate a new tree to be added onto the tip of the original tree, based on the function arguments specified.
    new.tree <- sim.bdtree(b=birth, d=death, stop="time", t=steps)
    # Rename tip labels of new tree to be prefixed by tip label on the original tree, where new tree will be bound.
    # Get number of tips of new tree.
    new.tips <- Ntip(new.tree)
    # For each tip in new tree,
    for (j in 1:new.tips){
      # get the current tip, and reformat its name.
      cur.tip <- new.tree$tip.label[j]
      # The seq.err flag denotes whether to add a "p" or an "e" to the new tips
      # (This is meant to represent either populations or sequencing error being added to the tree, respectively)
      if(type=="seq.err"){
        cur.tip <- sub("s","e",cur.tip) 
      } else{
        cur.tip <- sub("s","p",cur.tip)
      }
      new.name <- paste0(orig.tip,"--",cur.tip)
      # Add that new name to the tip of the new tree.
      new.tree$tip.label[j] <- new.name
    }
    # Bind newly generated tree to current tip of original tree.
    tree <- bind.tree(tree,new.tree,where=1)
    # The where argument of this command always equals 1, because every time a new tree is appended to the original, 
    # the next "open" tip of the original tree (which is where the next new tree will be added) becomes tip number 1.
  }
  # Drop extinct branches of the tree, such that only extant leaves are shown
  tree <- drop.extinct(tree)
  return(tree)
}

# -----------------------------------------
# %%% MPD VERSUS DELTA MAPPING FUNCTION %%%
# -----------------------------------------
# Function for capturing mpd values over series of delta transformations
phy.d.transform <- function(phylo,d,abundance.matrix,title){
  # phylo--the phylogenetic tree to be transformed
  # d--the vector of delta values to utilize for the branch length transformation
  # abundance.matrix--species-sites abundance matrix, used to calculate mpd values
  # title--title of plot generated (plotting mpd values against increasing delta values)
  # -----------------------------------------------------------------------------
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
    c.data <- comparative.comm(s.phylo, abundance.matrix)
    # Calculate mean pairwise distance
    mpd.values <- .mpd(c.data, abundance.weighted=TRUE) # This means that abundance data in simulation will effect calculation of mpd...Does this makes sense, given what we're trying to do?
    # Store calculated mpd values within a matrix
    mpd.mat[,i] <- mpd.values
  }
  # PLOTTING--For each 'site' or simulation, plot the change in mpd versus the increasing values of delta
  # Specifying the range values, removing any NAs
  #ymin <- min(mpd.mat, na.rm=T); ymax <- max(mpd.mat,na.rm=T)
  # Plotting the first matrix row (i.e. Site1 values) 
  #plot((mpd.mat[1,1:length(d)]) ~ d, xlab="delta", ylab="MPD Values", ylim=c(ymin,ymax), pch=20,main=title)
  #lines(d, mpd.mat[1,])
  # Below loop iterates through length of the matrix, adding connected points onto the plot
  #for(i in 2:nsim){
    #points((mpd.mat[i,1:length(d)]) ~ d,col=i, pch=20)
    # Capture mpd values for current row of matrix, and connect points
    #y <- (mpd.mat[i,])
    #lines(d,y,col=i)
  #}
  return(mpd.mat)
}

# ---------------------------------------
# %%% GENERATE VECTOR OF DELTA VALUES %%%
# ---------------------------------------
deltas <- seq(0.1,3,by=0.1)

# ---------------------------
# %%% SIMULATION FUNCTION %%%
# ---------------------------

# ABUNDANCE MAPPING FUNCTION
abundance.mapping <- function(oldColumnNames, newColumnNames, rowNames, donor.comm, recip.comm){
  # For each column in the donor matrix,
  for(i in 1:length(oldColumnNames)){
    # and for each site in the donor matrix,
    for(j in 1:length(rowNames)){
      # Get the columns in the recipient matrix that are derived from the "current" column of the donor matrix
      current.species <- grep(oldColumnNames[i],newColumnNames)
      # If no existing columns are derived from the "current" column in the donor 
      # (i.e. that species went extinct), then there are no values to copy, and we skip over that column
      if(length(current.species) == 0){
        next
      } else {
        # If there's only one remaining column derived from the "current" column in the donor, 
        # that column will receive the species' abundance values
        # (This if statement is necessary because of the unintuitive behavior of `sample` for vectors of length 1)
        if (length(current.species) == 1){
          lucky.column <- current.species
        } else {
          # If there are multiple columns derived from the "current" column, randomly select one of them
          lucky.column <- sample(x=current.species, size=1)
        }
        # The abundance values in the chosen ("lucky") column get the donor abundance values
        recip.comm[j,lucky.column] <- donor.comm[j,i]
        # Remove the "lucky.column" from the vector of populations of the current species
        # (This is to determine which remaining populations get abundances of 0)
        current.species <- current.species[!current.species %in% lucky.column]
        # The remaining columns get abundance values of 0
        recip.comm[j,current.species] <- 0
      }
    }
  }
  return(recip.comm)
}

# %%% FUNCTION FOR TESTING TRANSFORMATIONS OVER COMMUNITY SIMULATIONS %%%
SimulateCommnunity <- function(comm.size,comm.spp,comm.timesteps,comm.migrate,comm.env,comm.abund,comm.stoch,comm.speciate,intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
  # SIMULATE AND CLEANUP COMMUNITY
  # Simulate community (using pez function sim.meta.phy.comm)
  DemoCom <- sim.meta.phy.comm(size=comm.size, n.spp=comm.spp, timesteps=comm.timesteps, p.migrate=comm.migrate, env.lam=comm.env, abund.lam=comm.abund, stoch.lam=comm.stoch, p.speciate=comm.speciate)
  # %%% "CLEANUP" SIMULATED COMMUNITY, AND SETUP COMMUNITY WITH ADDED POPULATIONS %%%
  # Drop sites with one or less species in them from community matrix  (by removing rows with only one nonzero value in them)
  # We do this because our mpd calculation is abundance weighted, and mpd can only be calculated between different species
  sim.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
  # Drop extinct species from community matrix (by removing columns with only 0s in them)
  sim.comm <- sim.comm[,apply(sim.comm,2,function(x) !all(x==0))]
  # Get the names of remaining species in community (after removing extinctions and 'singletons')
  species.names <- colnames(sim.comm)
  # Trim phylogeny to only contain remaining species (using drop.tip)
  # (Tips are determined by looking for the tip.labels that are not in the species names)
  sim.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
  # (...need to account for scenarios in which all species are dropped (i.e. all species are extinct)?)
  # Get the site names
  sites <- rownames(sim.comm)
  
  # INTRASPECIFIC
  # Add intraspecific difference to simulated community phylogeny
  # (Extinct species are dropped from the tree within the add.branch function)
  intra.phy <- add.branch(sim.phy,birth=intra.birth,death=intra.death,steps=intra.steps,"pops")
  # Create an empty matrix for the intra.phy
  intra.comm <- matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label))
  population.names <- colnames(intra.comm)
  # Map abundance values from original community matrix to columns of intra.comm 
  intra.comm <- abundance.mapping(species.names, population.names, sites, sim.comm, intra.comm)
  
  # SEQ. ERROR
  # Add random error to intra.tree phylogeny
  seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,steps=seq.steps,"seq.err")
  # Create an empty matrix for the seq.phy
  seq.comm <- matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label))
  individual.names <- colnames(seq.comm)
  # Map abundance values from original community matrix to columns of seq.comm
  seq.comm <- abundance.mapping(species.names, individual.names, sites, sim.comm, seq.comm)
  
  # CAPTURE MPD VALUES OVER TRANSFORMATIONS FOR EACH COMMUNITY TYPE
  # Plotting phylogenetic diversity versus different delta transformations for original simulated community
  sim.test <- phy.d.transform(sim.phy, deltas, sim.comm,"Original community")
  # Plotting phylogenetic diversity versus different delta transformations for community with appended populations
  intra.test <- phy.d.transform(intra.phy, deltas, intra.comm,"Community with intraspecific differences added")
  # Seq. error
  seq.test <- phy.d.transform(seq.phy, deltas, seq.comm,"Community with seq. error differences added")
  
  # PACKAGE SIMULATION DATA
  # Export a list containing all simulation data: phylogenies, abundances, and mpd.mat for each community "type" (original, intra, and seq)
  community.abundances <- list(orig.community=sim.comm,intra.community=intra.comm,seq_err.community=seq.comm)
  community.phylogenies <- list(orig.phylo=sim.phy,intra.phylo=intra.phy,seq_err.phylo=seq.phy)
  community.transforms <- list(orig.transform=sim.test,intra.transform=intra.test,seq_err.transform=seq.test)
  simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms)
  return(simulation.data)
}

# -----------------------------------------------------

#sim.data <- SimulateCommnunity(comm.size=10, comm.spp=1,comm.timesteps=100,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=0.1,intra.death=0.5,intra.steps=5,seq.birth=0.1,seq.death=0.5,seq.steps=5)
#str(sim.data)

# Varying intraspecific and sequencing error rates only
params <- data.frame(expand.grid(comm.size=10,comm.spp=1,comm.timesteps=100,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),intra.steps=1:5,seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1),seq.steps=1:5))
#nrow(params)

#params <- data.frame(expand.grid(comm.size=10,comm.spp=1,comm.timesteps=100,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=0.1,intra.death=0.5,intra.steps=1:5,seq.birth=0.1,seq.death=0.5,seq.steps=1:5))
nrow(params)

sim.Results <- mcMap(function(i) SimulateCommnunity(params$comm.size[i],params$comm.spp[i],params$comm.timesteps[i],params$comm.migrate[i],params$comm.env[i],params$comm.abund[i],params$comm.stoch[i],params$comm.speciate[i],params$intra.birth[i],params$intra.death[i],params$intra.steps[i],params$seq.birth[i],params$seq.death[i],params$seq.steps[i]),1:nrow(params),mc.cores=12)

#write.csv(sim.Results,"~/OTU_Phylogeny/simResults.csv") 
save.image()
# After mcMap, use sapply with the measureShifts function on the mcMap output to do analysis.