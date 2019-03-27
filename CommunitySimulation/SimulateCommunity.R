library(pez)
library(ape)
library(geiger)

# %%% SIMULATE A COMMUNITY %%%
#DemoCom <- sim.meta.phy.comm(size=10, n.spp=1,timesteps=100, p.migrate = 0.02, env.lam = 10, abund.lam = 4, stoch.lam = 1, p.speciate = 0.06)
DemoCom <- sim.meta.phy.comm(size=10, n.spp=1,timesteps=50, p.migrate = 0.02, env.lam = 10, abund.lam = 4, stoch.lam = 1, p.speciate = 0.06)
plot(DemoCom$phy)

# %%% "CLEANUP" SIMULATED COMMUNITY, AND SETUP COMMUNITY WITH ADDED POPULATIONS %%%
# Drop sites with one or less species in them from community matrix  (by removing rows with only one nonzero value in them)
# We do this because our mpd calculation is abundance weighted, and mpd can only be calculated between different species
# (...does it make sense for our mpd to be abundance weighted? Think about what this means.)
sim.comm <- DemoCom$comm[apply(DemoCom$comm, 1, function(x) sum(x!=0) > 1),,drop=F]
# Drop extinct species from community matrix (by removing columns with only 0s in them)
sim.comm <- sim.comm[,apply(sim.comm,2,function(x) !all(x==0))]
# Get the names of remaining species in community (after removing extinctions and 'singletons')
species.names <- colnames(sim.comm)
sim.comm
# Trim phylogeny to only contain remaining species (using drop.tip)
# (Tips are determined by looking for the tip.labels that are not in the species names)
sim.phy <- drop.tip(DemoCom$phy, tip = DemoCom$phy$tip.label[!DemoCom$phy$tip.label %in% species.names])
plot(sim.phy)
# (...need to account for scenarios in which all species are dropped (i.e. all species are extinct)?)
# Get the site names
sites <- rownames(sim.comm)
# Add intraspecific difference to simulated community phylogeny
# (Extinct species are dropped from the tree within the add.branch function)
intra.phy <- add.branch(sim.phy,birth=0.1,death=0.1,steps=3,"pops")
plot(intra.phy)
# Create an empty matrix for the intra.tree phylogeny
intra.comm <- matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label))
population.names <- colnames(intra.comm)

# %%% MAP ABUNDANCES VALUES FROM COMMUNITY MATRIX TO RELEVANT POPULATIONS OF NEW MATRIX %%%
# For each column in the original community matrix,
for(i in 1:length(species.names)){
  # and for each site in the original community matrix,
  for(j in 1:length(sites)){
    # Get the columns in the intra.comm matrix that are derived from the "current" column of the original
    current.species <- grep(species.names[i],population.names)
    # If no existing populations are derived from the "current" column in the original 
    # (i.e. that species went extinct), then there are no values to copy, and we skip over that column
    if(length(current.species) == 0){
      next
    } else {
        # If there's only one remaining population derived from the "current" column in the original, 
        # that population will receive the species' abundance values
        # (This if statement is necessary because of the unintuitive behavior of `sample` for vectors of length 1)
        if (length(current.species) == 1){
          lucky.population <- current.species
        } else {
        # If there are multiple populations derived from the "current" column, randomly select one of them
        lucky.population <- sample(x=current.species, size=1)
        }
        # The abundance values in the chosen ("lucky") population get the original abundance values
        intra.comm[j,lucky.population] <- sim.comm[j,i]
        # Remove the "lucky.population" from the vector of populations of the current species
        # (This is to determine which remaining populations get abundances of 0)
        current.species <- current.species[!current.species %in% lucky.population]
        # The remaining populations get abundance values of 0
        intra.comm[j,current.species] <- 0
    }
  }
}
sim.comm
intra.comm

# %%% CAPTURE MPD VALUES FOR BOTH COMMUNITIES OVER DELTA TRANSFORMATIONS %%%
# Plotting phylogenetic diversity versus different delta transformations for original simulated community
sim.test <- phy.d.transform(sim.phy, deltas, sim.comm,"Original community")
sim.test

# Plotting phylogenetic diversity versus different delta transformations for community with appended populations
intra.test <- phy.d.transform(intra.phy, deltas, intra.comm,"Community with intraspecific differences added")
intra.test
