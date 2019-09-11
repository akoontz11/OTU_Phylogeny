# %%% MEASUREMENT OF MPD SHIFTS %%%

# 1. Average absolute change in site "ranking" (least diverse -> most diverse) versus changes in delta value
# 2. Measure Pearson correlation coefficient between the average change in site rankings and the delta values

measureRankingShifts <- function(orig.transform, deltas){
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
  # 2. --MEASURE PEARSON CORRELATION COEFFICIENT BETWEEN RANKING SHIFTS AND DELTA VALUES--
  r.value <- cor(meanRankChanges,deltaDifferences)
  # Return two measures and vector of differences in delta values
  measures <- list(mpdRankingChanges = meanRankChanges, r = r.value, deltaDifferences = deltaDifferences)
  return(measures)
}

# Demonstration with simulated species-level transforms 
orig.site <- measureRankingShifts(sim.Results[[4]]$transforms$orig.transform,deltas)
str(orig.site)
# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(orig.site$mpdRankingChanges[11:29] ~ orig.site$deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2),main="Original phylogenies: changes in diversity versus changes in delta values",xlab="Differences between delta values",ylab="mpdRankChanges");lines(orig.site$deltaDifferences[11:29], orig.site$mpdRankingChanges[11:29],col=1);points(orig.site$mpdRankingChanges[1:10] ~ orig.site$deltaDifferences[1:10]);lines(orig.site$deltaDifferences[1:10], orig.site$mpdRankingChanges[1:10],col=2)

# Demonstration with simulated population-level transforms 
intra.site <- measureRankingShifts(sim.Results[[4]]$transforms$intra.transform,deltas)
str(intra.site)
# Plot changes in mpd ranking of sites versus changes in delta values, 
plot(intra.site$mpdRankingChanges[11:29] ~ intra.site$deltaDifferences[11:29],ylim=c(0,17),xlim=c(0,2.2),main="Intra phylogenies: changes in diversity versus changes in delta values",xlab="Differences between delta values",ylab="mpdRankChanges");lines(intra.site$deltaDifferences[11:29], intra.site$mpdRankingChanges[11:29],col=1);points(intra.site$mpdRankingChanges[1:10] ~ intra.site$deltaDifferences[1:10]);lines(intra.site$deltaDifferences[1:10], intra.site$mpdRankingChanges[1:10],col=2)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Average differences in mpd values at adjacent delta transform values versus changes in delta value
# 2. Measure Pearson correlation coefficient between the average change in site rankings and the delta values

measureShifts <- function(orig.transform, deltas){
  # (This function assumes that the deltas vector is structured
  # such that delta = 1.0 is the 10th value in the vector, and values increment by 0.1)
  # 1. --BUILD A VECTOR OF AVERGAGE CHANGES IN DIVERSITY VALUES BETWEEN ADJACENT DELTA VALUES--
  # Generate a matrix to capture changes in diversity values between adjacent delta values
  div.changes <- matrix(NA, nrow=nrow(orig.transform), ncol=length(deltas)-1)
  # Generate a numeric vector of means of mpd shifts (between adjacent delta values)
  mean.changes <- vector(length=length(deltas)-1)
  # For each value in the vector of delta values
  for(i in 1:length(deltas)-1){
    # calculate the differences between mpd values of the next delta value and the current delta value
    div.changes[,i] <- orig.transform[,i+1] - orig.transform[,i]
    # Calculate the means of those differences, and pass that value into a vector
    mean.changes[i] <- mean(div.changes[,i])
  }
  # 2. --MEASURE PEARSON CORRELATION COEFFICIENT BETWEEN MPD VALUES AND DELTA VALUES--
  # Create a vector of equal length to the vector of means (of mpd changes)
  small.deltas <- deltas[1:29]
  # Calculate r.value
  r.value <- cor(mean.changes, small.deltas)
  # Return two measures and vector of differences in delta values
  measures <- list(meansOfDiversityChanges = mean.changes, r = r.value, small.deltas = small.deltas)
  return(measures)
}

test.site <- measureShifts(sim.Results[[4]]$transforms$orig.transform, deltas)
test.site$meansOfDiversityChanges

test.site <- measureShifts(sim.data$transforms$orig.transform, deltas)
test.site$meansOfDiversityChanges

plot(test.site$meansOfDiversityChanges ~ test.site$small.deltas)

lm(MPD ~ delta + MPD_true + n_spp, data=data)
# This requires mpd measurements, delta values, mpd values at delta = 1.0, and the number of species present in the community
# Is this for all sites, that is, MPD values for all sites, species for all sites? How is this structured? 

# Afterwards, how to apply this to all simulation settings...