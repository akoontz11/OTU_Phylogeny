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

# Ultrametric-only scenarios--------------------------------------------------
# 
# # Varying intraspecific/sequencing error birth/death rates, and number of community species
# params <- data.frame(expand.grid(comm.spp=c(5,10,15),intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1)))
# 
# # # Subsetting to only include instances in which birth >= death parameters 
# params <- subset(params, intra.birth >= intra.death)
# params <- subset(params, seq.birth >= seq.death)
# 
# # Running simulation on 12 cores
# sim.Results <- mcMap(function(i) SimulateCommnunity(comm.size=10,params$comm.spp[i],comm.timesteps=40,
#                                                    comm.migrate=0.02,comm.env=10,comm.abund=4,
#                                                    comm.stoch=1,comm.speciate=0.06,params$intra.birth[i],
#                                                    params$intra.death[i],intra.steps=3,params$seq.birth[i],
#                                                    params$seq.death[i],seq.steps=3),1:nrow(params),mc.cores=12)
# 
# # Saving results
# save.image("simResults/simResults.ULTRA_ONLY.RData")

# (Non-)ultrametric scenarios--------------------------------------------------
# 
# # Parameter set for ultrametric phylogenies (no positive death rates)
# ultra.params <- data.frame(expand.grid(comm.spp=c(5,10,15),
#                                 comm.birth=seq(0.1,0.5,0.1),
#                                 intra.birth=seq(0.1,0.5,0.1),seq.birth=seq(0.1,0.5,0.1)))
# 
# # Parameter set for non-ultrametric phylogenies (death rates)
# params <- data.frame(expand.grid(comm.spp=c(5,10,15),
#                                  comm.birth=seq(0.1,0.5,0.1),comm.death=c(0.1,0.2),
#                                  intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),
#                                  seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1)))
# # Subsetting to only include instances in which birth >= death parameters
# params <- subset(params, comm.birth >= comm.death)
# params <- subset(params, intra.birth >= intra.death)
# params <- subset(params, seq.birth >= seq.death)
# 
# 
# # Simulations for ultrametric scenarios
# sim.ultraResults <- mcMap(function(i) SimulateCommnunity(ultra.params$comm.spp[i],comm.size=10,ultra.params$comm.birth[i],
#                                                     comm.death=0,comm.env=1, comm.abund=1,
#                                                     ultra.params$intra.birth[i],intra.death=0,intra.steps=3,
#                                                     ultra.params$seq.birth[i],seq.death=0,seq.steps=3),
#                                                     1:nrow(ultra.params),mc.cores=12)
# 
# # Simulations for non-ultrametric scenarios
# sim.Results <- mcMap(function(i) SimulateCommnunity(params$comm.spp[i],comm.size=10,params$comm.birth[i],
#                                                     params$comm.death[i],comm.env=1, comm.abund=1,
#                                                     params$intra.birth[i],params$intra.death[i],intra.steps=3,
#                                                     params$seq.birth[i],params$seq.death[i],seq.steps=3),
#                                                     1:nrow(params),mc.cores=12)
# 
# 
# save.image("simResults/simResults.RData")

# (Non-)ultrametric scenarios, SESmpd--------------------------------------------------

# Non-ultrametric
# Varying inter/intra/seq birth/death rates
params <- data.frame(expand.grid(comm.spp=seq(50,100,200), 
                                 intra.birth=seq(0.1,0.3,0.1), intra.death=seq(0.1,0.2,0.1),
                                 seq.birth=seq(0.1,0.3,0.1), seq.death=seq(0.1,0.2,0.1)))

# Subsetting to only include instances in which birth > death parameters
params <- subset(params, intra.birth > intra.death)
params <- subset(params, seq.birth > seq.death)

# Running simulation
sim.Results <- mcMap(function(i) SimulateCommunity(params$comm.spp[i], comm.size=50,
                                                   params$intra.birth[i], params$intra.death[i],
                                                   params$seq.birth[i], params$seq.death[i],
                                                   tree.ratio=c(25,1,0.5)),
                     1:nrow(params), mc.cores=20)

# Ultrametric
# Varying inter/intra/seq birth rates
ultra.params <- data.frame(expand.grid(comm.spp=seq(50,100,200),
                                       intra.birth=seq(0.1,0.3,0.1),
                                       seq.birth=seq(0.1,0.3,0.1)))

# Running simulation
sim.ultra.Results <- mcMap(function(i) SimulateCommunity(ultra.params$comm.spp[i],comm.size=50,
                                                         inter.birth=1,inter.death=0,
                                                         ultra.params$intra.birth[i],intra.death=0,
                                                         ultra.params$seq.birth[i],seq.death=0,
                                                         tree.ratio=c(25,1,0.5)),
                           1:nrow(ultra.params), mc.cores=20)

# Saving results
save.image("simResults/simResults_20210511.RData")
