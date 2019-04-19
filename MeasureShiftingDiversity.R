# %%% MEASUREMENT OF MPD SHIFTS %%%

# 1. Change in site "ranking" (of least diverse -> most diverse, say) versus changes in delta value
# 2. Measure Pearson correlation coefficient between the mpd values and the delta values

# Function needs to be able to accomodate both of the above metrics for all 3 levels of 
# species diversity (i.e. original, with populations, and with sequencing error)

measureShifts <- function(orig.transform, deltas){
  # (This function assumes the structure of the deltas vector, 
  # such that delta = 1.0 is the 10th value in the vector, and values increment by 0.1)
  # 1. --BUILD VECTOR OF CHANGES IN STIE RANKING BY MPD--
  # Generate a numeric vector of site rankings (at delta=1.0)
  baseline.rank <- as.numeric(rank(sort(orig.transform[,10],decreasing = F)))
  # Generate a vector of names of these site rankings (at delta=1.0)
  baseline.names <- names(sort(orig.transform[,10],decreasing = F)) 
  # Generate a vector to capture mean changes between delta values
  meanRankChanges <- vector(length = length(deltas)-1)
  # Generate a vector to capture delta value differences
  deltaDifferences <- vector(length = length(deltas)-1)
  # For each value in the vector of delta values
  for(i in 1:length(deltas)-1){
    # Get the names of sites sorted from lowest to highest
    comparison <- names(sort(orig.transform[,i],decreasing = F))
    # Generate a numeric vector corresponding to how site rankings have changed from "baseline" delta value
    # (delta=1.0) to "current" delta value (delta[i])
    shifts <- match(baseline.names, comparison)
    # Calculate changes between two different rankings
    changes <- abs(shifts - baseline.rank)
    # Calculate mean changes, and place into vector
    meanRankChanges[i] <- mean(changes)
    # Calculate difference between current delta value and baseline delta value (delta=1.0)
    deltaDifferences[i] <- abs(deltas[i]-deltas[10])
  }
  # 2. --MEASURE PEARSON CORRELATION COEFFICIENT BETWEEN MPD AND DELTA VALUES--
  r.value <- cor(meanRankChanges,deltaDifferences)
  measures <- list(mpdRankingChanges = meanRankChanges, r =r.value)
  return(measures)
}
test.site <- measureShifts(sim.data$transforms$orig.transform,deltas)

# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(test.site$mpdRankingChanges[11:29] ~ deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2));lines(deltaDifferences[11:29], test.site$mpdRankingChanges[11:29],col=1);points(test.site$mpdRankingChanges[1:10] ~ deltaDifferences[1:10]);lines(deltaDifferences[1:10], test.site$mpdRankingChanges[1:10],col=2)
