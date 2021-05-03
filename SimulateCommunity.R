# SIMULATE COMMUNITIES AND THEIR PHYLOGENIES OVER DELTA TRANSFORMATIONS

library(pez)
library(ape)
library(geiger)

# %%% Function mapping abundance values from one community matrix to a newer, expanded matrix %%%----
abundance.mapping <- function(oldColumnNames, newColumnNames, rowNames, donor.comm, recip.comm){
  # If `donor.com` argument is empty, column names will be NULL; return the recip.comm argument as NULL
  if(is.null(oldColumnNames)){
    recip.comm <- NULL
    return(recip.comm)
  }else{
    # For each column in the donor matrix, and each site, get recipient matrix column derived from "current" donor matrix column
    for(i in 1:length(oldColumnNames)){
      for(j in 1:length(rowNames)){
        current.species <- grep(oldColumnNames[i],newColumnNames)
        # If no existing columns are derived from the "current" column in the donor (i.e. extinction), skip over that column
        if(length(current.species) == 0){
          next
        } else {
          # If there's only one derived column, that column receives abundance values
          # (This if statement is necessary because of the unintuitive behavior of `sample` for vectors of length 1...)
          if (length(current.species) == 1){
            lucky.column <- current.species
          } else {
            # If there are multiple columns derived from the "current" column, randomly select one of them
            lucky.column <- sample(x=current.species, size=1)
          }
          # The abundance values in the chosen ("lucky") column get the donor abundance values
          recip.comm[j,lucky.column] <- donor.comm[j,i]
          # Remove the "lucky.column" from vector of populations. Remaining columns get abundances of 0
          current.species <- current.species[!current.species %in% lucky.column]
          recip.comm[j,current.species] <- 0
        }
      }
    }
    return(recip.comm)
  }
}

# %%% Function generating original (interspecific) phylogenies and communities %%%----
sim.comm <- function(nspp=20, nsite=50, birth=1, death=0, tree.ratio=5){
  # Setup tree and environmental gradient
  tree <- sim.bdtree(n=nspp, b=birth, d=death)
  nspp <- length(tree$tip.label)
  intercept <- rnorm(nspp, sd=5)
  slope <- rnorm(nspp, sd=5)
  env <- (seq_len(nsite)-1) / (sd(seq_len(nsite)-1))
  
  # Rescale interspecific tree according to tree.ratio parameter, to control for ratio of inter/intra/seq variation
  tree$edge.length <- tree$edge.length*(tree.ratio/max(branching.times(tree)))

  # Create community values, based on nsite and env.str arguments
  # This is over-involved; a simpler generation of comm would be more suitable
  comm <- outer(intercept, rep(1,nsite)) + outer(slope, env)
  if(any(comm < 0)){
    comm[comm < 0] <- 0
  }
  comm <- rpois(prod(dim(comm)),as.numeric(comm))

  # Format and return
  comm <- t(matrix(comm, nrow=nspp, ncol=nsite, dimnames=list(tree$tip.label, seq_len(nsite)-1)))
  data <- comparative.comm(tree, comm)
  return(data)
}

# %%% Simulating communities: sim.comm, SESmpd %%%----
SimulateCommunity <- function(comm.spp,comm.size,inter.birth,inter.death,
                              intra.birth,intra.death,seq.birth,seq.death,
                              tree.ratio){
  # Simulate community (using sim.comm function); loop to account for NULL results
  for(i in 1:20){
    # Generate NULL for a DemoCom object that is in error
    DemoCom <- tryCatch(sim.comm(nspp=comm.spp,nsite=comm.size,birth=inter.birth,death=inter.death,tree.ratio=tree.ratio),error=function(cond) {return(NULL)})
    # If DemoCom is NULL, reattempt 10 times
    if(is.null(DemoCom)){
      max.iter <-10
      for(j in 1:max.iter){
        DemoCom <- tryCatch(sim.comm(nspp=comm.spp,nsite=comm.size,birth=inter.birth,death=inter.death,tree.ratio=tree.ratio),error=function(cond) {return(NULL)})
        if(!is.null(DemoCom)){
          break
        } else {
          if(j==max.iter){
            break
          }
        }
      }
    }
    # If after 10 times DemoCom is still NULL, assign NULL to other variables
    if(is.null(DemoCom)){
      orig.phy <- NULL ; orig.comm <- NULL ; species.names <- NULL ; sites <- NULL
      # Otherwise, cleanup community and assign populations
    } else {
      # Drop sites with one or less species in them (because mpd needs to measure communities)
      orig.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
      # Drop extinct species from community matrix, and get names of remaining species
      orig.comm <- orig.comm[,apply(orig.comm,2,function(x) !all(x==0))]
      species.names <- colnames(orig.comm)
      sites <- rownames(orig.comm)
      # Trim phylogeny to only contain remaining species (using drop.tip), and get site names
      orig.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
    }
    
    # If, after trimming phylogeny and cleaning community, phylogeny is NULL, assign NULL to community as well
    if(is.null(orig.phy)){
      orig.comm <- NULL ; species.names <- NULL ; sites <- NULL
    } else {
      # If orig.phy has two or less species, SESmpd will return NULL. Repeat, with more species
      if(length(species.names) >= 2){
        next
      } else {
        # Otherwise, break out of larger loop
        break
      }
    }
  }
  # If orig.phy length no longer conforms to specified tree.ratio, scale it
  if(max(branching.times(orig.phy)) != tree.ratio){
    orig.phy$edge.length <- orig.phy$edge.length*(tree.ratio/max(branching.times(orig.phy)))
  }
  
  # Add intraspecific differences----
  intra.phy <- add.branch(orig.phy,birth=intra.birth,death=intra.death,"pops")
  # Create an empty matrix for the intra community
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
  population.names <- colnames(intra.comm)
  # Map abundance values from original community matrix to columns of intra.comm
  intra.comm <- abundance.mapping(species.names, population.names, sites, orig.comm, intra.comm)
  
  # Add sequencing error----
  seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,"seq.err")
  # Create an empty matrix for the seq. error community
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
  individual.names <- colnames(seq.comm)
  # Map abundance values from original community matrix to columns of seq.comm
  seq.comm <- abundance.mapping(species.names, individual.names, sites, orig.comm, seq.comm)
  
  # Capture SESmpd values over transformations----
  # Original community
  orig.transform <- phy.d.transform(orig.phy, orig.comm, deltas)
  # Community with populations (intraspecific error)
  intra.transform <- phy.d.transform(intra.phy, intra.comm, deltas)
  # Community with sequencing error
  seq.transform <- phy.d.transform(seq.phy, seq.comm, deltas)
  
  # Capture mpd of original phylogenies (for later comparison)----
  orig.SESmpd <- tryCatch(.ses.mpd(comparative.comm(orig.phy, orig.comm, force.root = 0), abundance.weighted=TRUE), error=function(cond) {return(NULL)})
  orig.SESmpd <- orig.SESmpd$mpd.obs.z
  # Change names to match format from phy.d.transform function
  if(!is.null(orig.SESmpd)){
    names(orig.SESmpd) <- paste("Site",1:nrow(orig.comm),sep="_")
  }
  # Export a list containing all simulation data, for each community "type" (original, intra, and seq)----
  # Abundances
  community.abundances <- list(orig.community=orig.comm,intra.community=intra.comm,seq.community=seq.comm)
  # Phylogenies
  community.phylogenies <- list(orig.phylo=orig.phy,intra.phylo=intra.phy,seq.phylo=seq.phy)
  # MPD matrices, from delta transforms
  community.transforms <- list(orig.transform=orig.transform,intra.transform=intra.transform,seq.transform=seq.transform)
  # Original MPD values
  original.diversityMetrics <- list(SESmpds=orig.SESmpd)
  # Return data
  simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms,values=original.diversityMetrics)
  return(simulation.data)
}
# # Ultrametric (no death)
# test <- SimulateCommunity(comm.spp=10, comm.size=10, inter.birth=0.5, inter.death=0, tree.ratio=25,
#                   intra.birth=0.5, intra.death=0,
#                   seq.birth=0.5, seq.death=0)
# 
# # Non-ultrametric interspecific tree
# SimulateCommunity(comm.spp=10, comm.size=10, inter.birth=0.5, inter.death=0.2, tree.ratio=25,
#                   intra.birth=0.5, intra.death=0,
#                   seq.birth=0.5, seq.death=0)
# 
# # Non-ultrametric intra/seq branches
# test <- SimulateCommunity(comm.spp=10, comm.size=10, inter.birth=0.5, inter.death=0, tree.ratio=25,
#                   intra.birth=0.5, intra.death=0.2,
#                   seq.birth=0.5, seq.death=0.2)
# 
# # Non-ultrametric inter/intra/seq (high death)
# SimulateCommunity(comm.spp=5, comm.size=10, inter.birth=0.5, inter.death=0.2, tree.ratio=25,
#                   intra.birth=0.1, intra.death=0.5,
#                   seq.birth=0.1, seq.death=0.5)

# # %%% Simulate communities using sim.comm %%%----
# SimulateCommnunity <- function(comm.spp,comm.size,comm.birth,comm.death,comm.env,comm.abund,
#                                intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
#   # Simulate community (using Will's sim.comm function)
#   # Generate NULL for a DemoCom object that is in error
#   DemoCom <- tryCatch(sim.comm(nspp=comm.spp,nsite=comm.size,birth=comm.birth,death=comm.death,env.str=comm.env, min.lam=comm.abund),error=function(cond) {return(NULL)})
#   # If DemoCom is NULL, reattempt 100 times
#   if(is.null(DemoCom)){
#     max.iter <-100
#     for(i in 1:max.iter){
#       DemoCom <- tryCatch(sim.comm(nspp=comm.spp,nsite=comm.size,birth=comm.birth,death=comm.death,env.str=comm.env, min.lam=comm.abund),error=function(cond) {return(NULL)})
#       if(!is.null(DemoCom)){
#         break
#       } else {
#         if(i==max.iter){
#           break
#         }
#       }
#     }
#   }
#   # If after 100 times DemoCom is still NULL, assign NULL to other variables
#   if(is.null(DemoCom)){
#     orig.phy <- NULL ; orig.comm <- NULL ; species.names <- NULL ; sites <- NULL
#     # Otherwise, cleanup community and assign populations
#   } else {
#     # Drop sites with one or less species in them (because mpd needs to measure communities)
#     orig.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
#     # Drop extinct species from community matrix, and get names of remaining species
#     orig.comm <- orig.comm[,apply(orig.comm,2,function(x) !all(x==0))]
#     # orig.comm <- DemoCom$comm
#     species.names <- colnames(orig.comm)
#     sites <- rownames(orig.comm)
#     # Trim phylogeny to only contain remaining species (using drop.tip), and get site names
#     orig.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
#     # orig.phy <- DemoCom$phy
#   }
# 
#   # Add intraspecific differences----
#   intra.phy <- add.branch(orig.phy,birth=intra.birth,death=intra.death,steps=intra.steps,"pops")
#   # Create an empty matrix for the intra.phy
#   # (tryCatch statement included to account for trees in which all species have gone extinct)
#   intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
#   population.names <- colnames(intra.comm)
#   # Map abundance values from original community matrix to columns of intra.comm
#   intra.comm <- abundance.mapping(species.names, population.names, sites, orig.comm, intra.comm)
# 
#   # Add sequencing error----
#   seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,steps=seq.steps,"seq.err")
#   # Create an empty matrix for the seq.phy
#   # (tryCatch statement included to account for trees in which all species have gone extinct)
#   seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
#   individual.names <- colnames(seq.comm)
#   # Map abundance values from original community matrix to columns of seq.comm
#   seq.comm <- abundance.mapping(species.names, individual.names, sites, orig.comm, seq.comm)
# 
#   # Capture mpd values over transformations----
#   # Original community
#   orig.transform <- phy.d.transform(orig.phy, orig.comm, deltas)
#   # Community with populations (intraspecific error)
#   intra.transform <- phy.d.transform(intra.phy, intra.comm, deltas)
#   # Community with sequencing error
#   seq.transform <- phy.d.transform(seq.phy, seq.comm, deltas)
# 
#   # Capture mpd of original phylogenies (for later comparison)----
#   orig.MPD <- tryCatch(.mpd(comparative.comm(orig.phy, orig.comm, force.root = 0), abundance.weighted=TRUE), error=function(cond) {return(NULL)})
#   # Change names to match format from phy.d.transform function
#   if(!is.null(orig.MPD)){
#     names(orig.MPD) <- paste("Site",1:nrow(orig.comm),sep="_")
#   }
# 
#   # Export a list containing all simulation data, for each community "type" (original, intra, and seq)----
#   # Abundances
#   community.abundances <- list(orig.community=orig.comm,intra.community=intra.comm,seq.community=seq.comm)
#   # Phylogenies
#   community.phylogenies <- list(orig.phylo=orig.phy,intra.phylo=intra.phy,seq.phylo=seq.phy)
#   # MPD matrices, from delta transforms
#   community.transforms <- list(orig.transform=orig.transform,intra.transform=intra.transform,seq.transform=seq.transform)
#   # Original MPD values
#   original.diversityMetrics <- list(MPDs=orig.MPD)
#   # Return data
#   simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms,values=original.diversityMetrics)
#   return(simulation.data)
# }

# # Ultrametric test
# SimulateCommnunity(comm.spp=10, comm.size=10, comm.birth=0.5, comm.death=0,
#                    comm.env=1, comm.abund=1,intra.birth=0.5,
#                    intra.death=0.1,intra.steps=3,seq.birth=0.5,
#                    seq.death=0.1,seq.steps=3)
# 
# # Non-ultrametric tests
# SimulateCommnunity(comm.spp=10, comm.size=10, comm.birth=0.5, comm.death=0.1,
#                            comm.env=1, comm.abund=1,intra.birth=0.5,
#                            intra.death=0.1,intra.steps=3,seq.birth=0.5,
#                            seq.death=0.1,seq.steps=3)
# # High death
# SimulateCommnunity(comm.spp=5, comm.size=10, comm.birth=0.5, comm.death=0.2,
#                           comm.env=1, comm.abund=1,intra.birth=0.1,
#                           intra.death=0.5,intra.steps=3,seq.birth=0.1,
#                           seq.death=0.5,seq.steps=3)
# plot(test$phylogenies$orig.phylo, show.tip.label = FALSE)
# plot(test$phylogenies$seq.phylo, show.tip.label = FALSE)
# test$values

# # %%% Simulate communities using sim.meta.phy.comm %%%----
# SimulateCommnunity <- function(comm.size,comm.spp,comm.timesteps,comm.migrate,comm.env,comm.abund,comm.stoch,comm.speciate,intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
# Simulate community (using pez function sim.meta.phy.comm)
#   # Generate NULL for a DemoCom object that is in error
#   DemoCom <- tryCatch(sim.meta.phy.comm(size=comm.size, n.spp=comm.spp, timesteps=comm.timesteps, p.migrate=comm.migrate, env.lam=comm.env, abund.lam=comm.abund, stoch.lam=comm.stoch, p.speciate=comm.speciate),error=function(cond) {return(NULL)})
#   # If DemoCom is NULL, then assign NULL to relevant phylogeny/community variables
#   if(is.null(DemoCom)){
#     orig.phy <- NULL ; orig.comm <- NULL ; species.names <- NULL ; sites <- NULL
#   }else{
#
#     # Cleanup community and assign populations
#     # Drop sites with one or less species in them (because mpd needs to measure communities)
#     sim.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
#     # Drop extinct species from community matrix, and get names of remaining species
#     orig.comm <- sim.comm[,apply(sim.comm,2,function(x) !all(x==0))]
#     species.names <- colnames(orig.comm)
#     # Trim phylogeny to only contain remaining species (using drop.tip), and get site names
#     orig.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
#     sites <- rownames(orig.comm)
#   }
#
# Add intraspecific differences
#   intra.phy <- add.branch(orig.phy,birth=intra.birth,death=intra.death,steps=intra.steps,"pops")
#   # Create an empty matrix for the intra.phy
#   # (tryCatch statement included to account for trees in which all species have gone extinct)
#   intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
#   population.names <- colnames(intra.comm)
#   # Map abundance values from original community matrix to columns of intra.comm
#   intra.comm <- abundance.mapping(species.names, population.names, sites, orig.comm, intra.comm)
#
# Add sequencing error
#   seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,steps=seq.steps,"seq.err")
#   # Create an empty matrix for the seq.phy
#   # (tryCatch statement included to account for trees in which all species have gone extinct)
#   seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
#   individual.names <- colnames(seq.comm)
#   # Map abundance values from original community matrix to columns of seq.comm
#   seq.comm <- abundance.mapping(species.names, individual.names, sites, sim.comm, seq.comm)
#
# Capture mpd values over transformations
#   # Original community
#   orig.test <- phy.d.transform(orig.phy, orig.comm, deltas)
#   # Community with populations (intraspecific error)
#   intra.test <- phy.d.transform(intra.phy, intra.comm, deltas)
#   # Community with sequencing error
#   seq.test <- phy.d.transform(seq.phy, seq.comm, deltas)
#
# Capture mpd of original phylogenies (for later comparison)
#   orig.MPD <- .mpd(comparative.comm(orig.phy, orig.comm, force.root = 0), abundance.weighted=TRUE)
#   # Change names to match format from phy.d.transform function
#   names(orig.MPD) <- paste("Site",1:nrow(orig.comm),sep="_")
#
# Export a list containing all simulation data, for each community "type" (original, intra, and seq)
#   # Abundances
#   community.abundances <- list(orig.community=orig.comm,intra.community=intra.comm,seq.community=seq.comm)
#   # Phylogenies
#   community.phylogenies <- list(orig.phylo=orig.phy,intra.phylo=intra.phy,seq.phylo=seq.phy)
#   # MPD matrices, from delta transforms
#   community.transforms <- list(orig.transform=orig.test,intra.transform=intra.test,seq.transform=seq.test)
#   # Original MPD values
#   original.diversityMetrics <- list(MPDs=orig.MPD)
#   # Return data
#   simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms,values=original.diversityMetrics)
#   return(simulation.data)
# }
# 
# 
# test <- SimulateCommnunity(comm.size=10,comm.spp=10,comm.timesteps=40,
#                            comm.migrate=0.02,comm.env=10,comm.abund=4,
#                            comm.stoch=1,comm.speciate=0.06,intra.birth=0.5,
#                            intra.death=0.1,intra.steps=3,seq.birth=0.5,
#                            seq.death=0.1,seq.steps=3)
