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
