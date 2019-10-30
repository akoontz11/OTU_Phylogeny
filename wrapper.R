library(pez)
library(ape)
library(geiger)
library(phytools)
library(parallel)

# Add.branch function
source("~/OTU_Phylogeny/Intraspecific_SeqError_Sim.R")
# MPD vs. delta function
source("~/OTU_Phylogeny/MPDvsDelta.R")
# Generate vector of delta values
deltas <- seq(0.1,3,by=0.1)
# Simulate community function
source("~/OTU_Phylogeny/SimulateCommunity.R")

# -----------------------------------------------------

# Varying intraspecific rates, sequencing error rates, and number of community species
params <- data.frame(expand.grid(comm.size=10,comm.spp=c(5,10,15),comm.timesteps=40,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=0.06,intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),intra.steps=1:5,seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1),seq.steps=1:5))
# Subsetting data.frame to only include instances in which birth >= death parameters (for both intra and seq)
params <- subset.data.frame(params, intra.birth >= intra.death)
params <- subset.data.frame(params, seq.birth >= seq.death)

# Running simulation on 12 cores
sim.Results <- mcMap(function(i) SimulateCommnunity(params$comm.size[i],params$comm.spp[i],params$comm.timesteps[i],params$comm.migrate[i],params$comm.env[i],params$comm.abund[i],params$comm.stoch[i],params$comm.speciate[i],params$intra.birth[i],params$intra.death[i],params$intra.steps[i],params$seq.birth[i],params$seq.death[i],params$seq.steps[i]),1:nrow(params),mc.cores=12)

# Saving results
save.image("simResults.RData")