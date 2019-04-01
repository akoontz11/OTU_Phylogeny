library(pez)
library(ape)
library(geiger)

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
  # (...does it make sense for our mpd to be abundance weighted? Think about what this means.)
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
  # Export a list of each of the matrices of mpd values at each site for each transform value
  transform.matrices <- list(sim.test,intra.test,seq.test)
  return(transform.matrices)
}

test.matrices <- SimulateCommnunity(comm.size=10, comm.spp=1,comm.timesteps=100,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=0.1,intra.death=0.1,intra.steps=3,seq.birth=0.2,seq.death=0.2,seq.steps=3)
str(test.matrices)
