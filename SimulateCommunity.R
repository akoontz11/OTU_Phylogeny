library(pez)
library(ape)
library(geiger)

# ABUNDANCE MAPPING FUNCTION
# Function for mapping abundance values from one community matrix to a newer, expanded matrix
abundance.mapping <- function(oldColumnNames, newColumnNames, rowNames, donor.comm, recip.comm){
  # If `donor.com` argument is empty, column names will be NULL; return the recip.comm argument as NULL
  if(is.null(oldColumnNames)){
    recip.comm <- NULL
    return(recip.comm)
  }else{
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
          # (This if statement is necessary because of the unintuitive behavior of `sample` for vectors of length 1...)
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
}

# %%% FUNCTION FOR TESTING TRANSFORMATIONS OVER COMMUNITY SIMULATIONS %%%
SimulateCommnunity <- function(comm.size,comm.spp,comm.timesteps,comm.migrate,comm.env,comm.abund,comm.stoch,comm.speciate,intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
  #browser()
  # SIMULATE AND CLEANUP COMMUNITY
  # Simulate community (using pez function sim.meta.phy.comm). Generate NULL for a DemoCom object that is in error
  DemoCom <- tryCatch(sim.meta.phy.comm(size=comm.size, n.spp=comm.spp, timesteps=comm.timesteps, p.migrate=comm.migrate, env.lam=comm.env, abund.lam=comm.abund, stoch.lam=comm.stoch, p.speciate=comm.speciate),error=function(cond) {return(NULL)})
  # If DemoCom is NULL, then assign NULL to relevant phylogeny/community variables
  if(is.null(DemoCom)){
    orig.phy <- NULL ; orig.comm <- NULL ; species.names <- NULL ; sites <- NULL 
  }else{
  # %%% "CLEANUP" SIMULATED COMMUNITY, AND SETUP COMMUNITY WITH POPULATIONS %%%
  # Drop sites with one or less species in them from community matrix  (by removing rows with only one nonzero value in them)
  # (We do this because our mpd calculation is abundance weighted, and mpd can only be calculated between different species)
  sim.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
  # Drop extinct species from community matrix (by removing columns with only 0s in them)
  orig.comm <- sim.comm[,apply(sim.comm,2,function(x) !all(x==0))]
  # Get the names of remaining species in community (after removing extinctions and 'singletons')
  species.names <- colnames(orig.comm)
  # Trim phylogeny to only contain remaining species (using drop.tip)
  # (Tips are determined by looking for the tip.labels that are not in the species names)
  orig.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
  # Get the site names
  sites <- rownames(orig.comm)
  }
  # INTRASPECIFIC
  # Add intraspecific difference to simulated community phylogeny
  intra.phy <- add.branch(orig.phy,birth=intra.birth,death=intra.death,steps=intra.steps,"pops")
  # Create an empty matrix for the intra.phy
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
  population.names <- colnames(intra.comm)
  # Map abundance values from original community matrix to columns of intra.comm 
  intra.comm <- abundance.mapping(species.names, population.names, sites, orig.comm, intra.comm)
  
  # SEQ. ERROR
  # Add random error to intra.tree phylogeny
  seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,steps=seq.steps,"seq.err")
  # Create an empty matrix for the seq.phy
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
  individual.names <- colnames(seq.comm)
  # Map abundance values from original community matrix to columns of seq.comm
  seq.comm <- abundance.mapping(species.names, individual.names, sites, sim.comm, seq.comm)
  
  # CAPTURE MPD VALUES OVER TRANSFORMATIONS FOR EACH COMMUNITY TYPE
  # Build mpd versus delta transformation values matrix for original simulated community
  orig.test <- phy.d.transform(orig.phy, orig.comm, deltas)
  # Build mpd versus delta transformation values matrix for community with appended populations
  intra.test <- phy.d.transform(intra.phy, intra.comm, deltas)
  # Build mpd versus delta transformation values matrix for community with sequencing error added on
  seq.test <- phy.d.transform(seq.phy, seq.comm, deltas)
  
  # CAPTURE PHYLOGENETIC DIVERSITY METRICS OF ORIGINALLY SIMULATED PHYLOGENIES (FOR LATER COMPARISON)
  # --Originally simulated phylogeny/community--
  # Mean pairwise distance
  orig.MPD <- .mpd(comparative.comm(orig.phy, orig.comm, force.root = 0), abundance.weighted=TRUE)
  # Ranking of sites by diversity
  names(orig.MPD) <- rownames(orig.comm)
  orig.ranking <- names(sort(orig.MPD,decreasing = F))
  # --Intra phylogeny/community--
  # Mean pairwise distance
  intra.MPD <- .mpd(comparative.comm(intra.phy, intra.comm, force.root = 0), abundance.weighted=TRUE)
  # Ranking of sites by diversity
  names(intra.MPD) <- rownames(intra.comm)
  intra.ranking <- names(sort(intra.MPD,decreasing = F))
  # --Seq phylogeny/community--
  # Mean pairwise distance
  seq.MPD <- .mpd(comparative.comm(seq.phy, seq.comm, force.root = 0), abundance.weighted=TRUE)
  # Ranking of sites by diversity
  names(seq.MPD) <- rownames(seq.comm)
  seq.ranking <- names(sort(seq.MPD,decreasing = F))
  # Combining metrics into lists
  MPDs <- list(orig=orig.MPD, intra=intra.MPD, seq=seq.MPD)
  Rankings <- list(orig=orig.ranking, intra=intra.ranking, seq=seq.ranking)
  
  # PACKAGE SIMULATION DATA
  # Export a list containing all simulation data, for each community "type" (original, intra, and seq)
  # Abundances
  community.abundances <- list(orig.community=orig.comm,intra.community=intra.comm,seq.community=seq.comm)
  # Phylogenies
  community.phylogenies <- list(orig.phylo=orig.phy,intra.phylo=intra.phy,seq.phylo=seq.phy)
  # MPD matrices, from delta transforms
  community.transforms <- list(orig.transform=orig.test,intra.transform=intra.test,seq.transform=seq.test)
  # Phylogenetic diversity metrics, of original (untransformed) communities/phylogenies
  community.diversityMetrics <- list(mpds=MPDs, rankings=Rankings)
  # Return data
  simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms,values=community.diversityMetrics)
  return(simulation.data)
}

sim.data <- SimulateCommnunity(comm.size=10, comm.spp=10,comm.timesteps=40,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=0.5,intra.death=0.1,intra.steps=1,seq.birth=0.5,seq.death=0.1,seq.steps=1)
str(sim.data)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMULATE AND CLEANUP COMMUNITY
# Simulate community (using pez function sim.meta.phy.comm)
# DemoCom <- tryCatch(sim.meta.phy.comm(size=5, n.spp=1, timesteps=100, p.migrate=0.02, env.lam=10, abund.lam=4, stoch.lam=1, p.speciate=0.01),error=function(cond) {return(NULL)})
# #DemoCom <- tryCatch(sim.meta.phy.comm(size=10, n.spp=1, timesteps=100, p.migrate=0.02, env.lam=10, abund.lam=4, stoch.lam=1, p.speciate=0.06),error=function(cond) {return(NULL)})
# # If DemoCom is NULL, then assign NULL to relevant variables
# if(is.null(DemoCom)){
#   sim.phy <- NULL ; sim.comm <- NULL ; species.names <- NULL ; sites <- NULL
# }else{
#   # %%% "CLEANUP" SIMULATED COMMUNITY, AND SETUP COMMUNITY WITH ADDED POPULATIONS %%%
#   # Drop sites with one or less species in them from community matrix  (by removing rows with only one nonzero value in them)
#   # (We do this because our mpd calculation is abundance weighted, and mpd can only be calculated between different species)
#   sim.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
#   # Drop extinct species from community matrix (by removing columns with only 0s in them)
#   sim.comm <- sim.comm[,apply(sim.comm,2,function(x) !all(x==0))]
#   # Get the names of remaining species in community (after removing extinctions and 'singletons')
#   species.names <- colnames(sim.comm)
#   # Trim phylogeny to only contain remaining species (using drop.tip)
#   # (Tips are determined by looking for the tip.labels that are not in the species names)
#   sim.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
#   # Get the site names
#   sites <- rownames(sim.comm)
# }
# # INTRASPECIFIC
# # Add intraspecific difference to simulated community phylogeny
# # (tryCatch statements included to account for trees in which all species have gone extinct)
# intra.phy <- add.branch(sim.phy,birth=0.4,death=0.2,steps=1,"pops")
# # Create an empty matrix for the intra.phy
# intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
# population.names <- colnames(intra.comm)
# # Map abundance values from original community matrix to columns of intra.comm
# intra.comm <- abundance.mapping(species.names, population.names, sites, sim.comm, intra.comm)
# 
# # SEQ. ERROR
# # Add random error to intra.tree phylogeny
# # (tryCatch statements included to account for trees in which all species have gone extinct)
# seq.phy <- add.branch(intra.phy,birth=0.1,death=0.1,steps=1,"seq.err")
# # Create an empty matrix for the seq.phy
# seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
# individual.names <- colnames(seq.comm)
# # Map abundance values from original community matrix to columns of seq.comm
# seq.comm <- abundance.mapping(species.names, individual.names, sites, sim.comm, seq.comm)
# 
# # CAPTURE MPD VALUES OVER TRANSFORMATIONS FOR EACH COMMUNITY TYPE
# # Plotting phylogenetic diversity versus different delta transformations for original simulated community
# sim.test <- phy.d.transform(sim.phy, deltas, sim.comm,"Original community")
# # Plotting phylogenetic diversity versus different delta transformations for community with appended populations
# intra.test <- phy.d.transform(intra.phy, deltas, intra.comm,"Community with intraspecific differences added")
# # Plotting phylogenetic diversity versus different delta transformations for community with sequencing error added on
# seq.test <- phy.d.transform(seq.phy, deltas, seq.comm,"Community with seq. error differences added")
# 
# # PACKAGE SIMULATION DATA
# # Export a list containing all simulation data: phylogenies, abundances, and mpd.mat for each community "type" (original, intra, and seq)
# community.abundances <- list(orig.community=sim.comm,intra.community=intra.comm,seq_err.community=seq.comm)
# community.phylogenies <- list(orig.phylo=sim.phy,intra.phylo=intra.phy,seq_err.phylo=seq.phy)
# community.transforms <- list(orig.transform=sim.test,intra.transform=intra.test,seq_err.transform=seq.test)
# simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms)
