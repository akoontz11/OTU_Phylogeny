library(parallel)

wrapper <- function(comm.size,comm.spp,comm.timesteps,comm.migrate,comm.env,comm.abund,comm.stoch,comm.speciate,intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
  # Simulate community, and add levels of population and sequencing error onto community and phylogeny.
  # Then, calculate mpd at each value of delta transformations and for all 3 levels of "genetic granularity"
  sim.data <- SimulateCommnunity(comm.size=comm.size, comm.spp=com.spp,comm.timesteps=comm.timesteps,comm.migrate=comm.migrate,comm.env=comm.env,comm.abund=comm.abund,comm.stoch=comm.stoch,comm.speciate=comm.speciate,intra.birth=intra.birth,intra.death=intra.death,intra.steps=intra.steps,seq.birth=seq.birth,seq.death=seq.death,seq.steps=seq.steps)
  # Measure change in mpd ranking and r value between mpd ranking change and delta transform value (different from 1)
  # --for original community
  baseline.community <- measureShifts(sim.data$transforms$orig.transform,deltas)
  # --for community with populations included
  population.community <- measureShifts(sim.data$transforms$intra.transform,deltas)
  #--for community with populations included and sequencing error added on
  sequenceErr.community <- measureShifts(sim.data$transforms$seq_err.transform,deltas)
  return() # --a list of 3 matrices and 3 vectors?
}

# measureShifts returns a list containing a vector (average changes in "mpd ranking" for each site in the community) 
# and a variable (r value between changes in mpd rank and changes in delta transformation values). 

# Since I want to capture the mpd ranking changes and the r value at all 3 levels of "species granularity" 
# (original community/phylogeny, intraspecific diffs, seq. error diffs), does this mean I need
# my wrapper to return a list of 3 matrices and 3 vectors?

# Building params through expand.grid
# The line below is huge--could we reduce this by removing some of these variables?
#params <- expand.grid(comm.size=1:5,comm.spp=1:5,comm.timesteps=100,comm.migrate=seq(0.01,0.05,0.01),comm.env=1:5,comm.abund=1:5,comm.stoch=1:5,comm.speciate=seq(0.01,0.05,0.01),intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),intra.steps=1:5,seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1),seq.steps=1:5)
# With comm.timesteps, comm.migrate, comm.env, comm.abund, and comm.stoch variables fixed:
params <- expand.grid(comm.size=1:5,comm.spp=1:5,comm.timesteps=100,comm.migrate=0.02,comm.env=10,comm.abund=4,comm.stoch=1,comm.speciate=seq(0.01,0.05,0.01),intra.birth=seq(0.1,0.5,0.1),intra.death=seq(0.1,0.5,0.1),intra.steps=1:5,seq.birth=seq(0.1,0.5,0.1),seq.death=seq(0.1,0.5,0.1),seq.steps=1:5)


output <- mcMap(SimulateCommnunity(),params,mc.cores=12)

library(parallel)
input <- list(1:10, letters[1:6], c(TRUE,FALSE,TRUE))
mclapply(input, length, mc.cores=2)
#...which is the same as...
mcMap(length, input, mc.cores=2)


# After mcMap, use sapply with the measureShifts function on the mcMap output to do analysis.