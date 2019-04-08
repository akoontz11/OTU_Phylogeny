# %%% MEASUREMENT OF MPD SHIFTS %%%

# 1. Change in site "ranking" (of least diverse -> most diverse, say) versus changes in delta value
  # a. Obtain mpd values for each corresponding delta value for an individual site
  # b. Get the average difference between mpd values for one delta value and mpd values for the next delta value
  # c. Plot those averages over differences between "current" and "baseline" delta values

# Generate a numeric vector of site rankings (at delta=1.0)
baseline.rank <- as.numeric(rank(sort(sim.data$transforms$orig.transform[,10],decreasing = F)))
# Generate a vector of names of these site rankings (at delta=1.0)
baseline.names <- names(sort(sim.data$transforms$orig.transform[,10],decreasing = F)) 
baseline.names
# Generate a vector to capture mean changes between delta values
meanRankChanges <- vector(length = length(deltas)-1)
# Generate a vector to capture delta value differences
deltaDifferences <- vector(length = length(deltas)-1)
# For each value in the vector of delta values
for(i in 1:length(deltas)-1){
  # Get the names of sites sorted from lowest to highest
  comparison <- names(sort(sim.data$transforms$orig.transform[,i],decreasing = F))
  # Generate a numeric vector corresponding to how site rankings have changed from "baseline" delta value
  # (delta=1.0) to "current" delta value (delta[i])
  shifts <- match(baseline.names, comparison)
  # Calculate changes between two different rankings
  changes <- abs(shifts - baseline.rank)
  # Calculate mean changes, and place into vector
  meanRankChanges[i] <- mean(changes)
  # Calculate difference between current delta value and baseline delta value (delta=1.0)
  deltaDifferences[i] <- (deltas[i]-deltas[10])
}

meanRankChanges
deltaDifferences

plot(meanRankChanges ~ deltaDifferences);lines(deltaDifferences, meanRankChanges)

# 2. Measure Pearson correlation coefficient between the mpd values and the delta values
my.model <- lm(meanRankChanges ~ deltaDifferences)
summary(my.model)

