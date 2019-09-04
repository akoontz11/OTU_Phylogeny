# %%% MEASUREMENT OF MPD SHIFTS %%%

# 1. Average absolute change in site "ranking" (least diverse -> most diverse) versus changes in delta value
# 2. Measure Pearson correlation coefficient between the mpd values and the delta values

# Function needs to be able to accomodate both of the above metrics for all 3 "levels" of 
# diversity (i.e. original, with populations, and with sequencing error)

measureShifts <- function(orig.transform, deltas){
  browser()
  # (This function assumes that the deltas vector is structured
  # such that delta = 1.0 is the 10th value in the vector, and values increment by 0.1)
  # 1. --BUILD VECTOR OF CHANGES IN STIE RANKING BY MPD--
  # Generate a numeric vector of site rankings (for delta=1.0)
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
  # Return two measures and vector of differences in delta values
  measures <- list(mpdRankingChanges = meanRankChanges, r = r.value, deltaDifferences = deltaDifferences)
  return(measures)
}
# Demonstration with simulated species-level transforms 
orig.site <- measureShifts(sim.data$transforms$orig.transform,deltas)
str(orig.site)
# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(orig.site$mpdRankingChanges[11:29] ~ orig.site$deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2),main="Original phylogenies: changes in diversity versus changes in delta values",xlab="Differences between delta values",ylab="mpdRankChanges");lines(orig.site$deltaDifferences[11:29], orig.site$mpdRankingChanges[11:29],col=1);points(orig.site$mpdRankingChanges[1:10] ~ orig.site$deltaDifferences[1:10]);lines(orig.site$deltaDifferences[1:10], orig.site$mpdRankingChanges[1:10],col=2)

# Demonstration with simulated population-level transforms
intra.site <- measureShifts(sim.data$transforms$intra.transform,deltas)
str(intra.site)
# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(intra.site$mpdRankingChanges[11:29] ~ intra.site$deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2),main="Intra- phylogenies: changes in diversity versus changes in delta values",xlab="Differences between delta values",ylab="mpdRankChanges");lines(intra.site$deltaDifferences[11:29], intra.site$mpdRankingChanges[11:29],col=1);points(intra.site$mpdRankingChanges[1:10] ~ intra.site$deltaDifferences[1:10]);lines(intra.site$deltaDifferences[1:10], intra.site$mpdRankingChanges[1:10],col=2)

# Demonstration with simulated seq-level transforms
seq.site <- measureShifts(sim.data$transforms$seq_err.transform,deltas)
str(seq.site)
# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(seq.site$mpdRankingChanges[11:29] ~ seq.site$deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2),main="Seq-error phylogenies: changes in diversity versus changes in delta values",xlab="Differences between delta values",ylab="mpdRankChanges");lines(seq.site$deltaDifferences[11:29], seq.site$mpdRankingChanges[11:29],col=1);points(seq.site$mpdRankingChanges[1:10] ~ seq.site$deltaDifferences[1:10]);lines(seq.site$deltaDifferences[1:10], seq.site$mpdRankingChanges[1:10],col=2)
