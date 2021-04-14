# ANALYSIS OF SIMULATION OUTPUT 

library(xtable)
library(viridis)

# %%% Worker functions %%%----
# Function to calculate correlations of MPD measurements between each delta value (x)
# and a "reference" delta value (d, taken prior to transformations)
.new.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}

# Function to calculate difference in site (MPD) rankings 
.ranking.diff <- function(x,d){
  # Generate a numeric vector of baseline site rankings (i.e. delta=1.0)
  baseline.rank <- rank(d)
  # Generate a vector to capture mean changes between rankings
  meanRankChanges <- vector(length = ncol(x))
  for (i in 1:ncol(x)){
    # Calculate the mean of the absolute changes between two different rankings, and place into vector
    meanRankChanges[i] <- mean(abs(baseline.rank - rank(x[,i])))
  }
  return(meanRankChanges)
}

# %%% Read in simulation data %%%----
load("simResults.RData")
# Generate backup data of simulation variables
backup <- sim.Results
b.params <- params

# Extract MPDs over delta transforms, from communities/phylogenies with intraspecific AND seq. error branches
sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
# Extract diversity metrics from untransformed phylogenies, for calculating response metrics
sim.MPDs <- lapply(sim.Results, function(x) x$values$MPDs)

# Build results matrix, from which linear model variables will be pulled
results <- params[rep(1:nrow(params), each=30),]
results$delta <- rep(deltas,length(sim.data))
# Create term for diversification (birth rate/death rate) for both intraspecific and seq.err.
results$intra.div <- results$intra.birth/results$intra.death
results$seq.div <- results$seq.birth/results$seq.death

# %%% MPD correlations between site and baseline %%%----
# Calculate MPD correlations on simulation data
t.correls <- mapply(.new.correls, sim.data, sim.MPDs)
# Match the correlations to the parameters
results$correl <- as.numeric(t.correls)
# Model effects on MPD correlations
s.model.correl <- lm(correl ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(s.model.correl)

# %%% Ranking differences between site and baseline %%%----
# Calculate mean number of site rank shifts
rank.shifts <- mapply(.ranking.diff, sim.data, sim.MPDs)
# Match the ranking differences to the parameters
results$rank.shifts <- as.numeric(rank.shifts)
# Model effects on site diversity rankings
s.model.rankShifts <- lm(rank.shifts ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(s.model.rankShifts)

# %%% Non-ultrametric runs %%%----
load("OTU_Phylogeny/simResults.20210412.RData")
backup <- sim.Results
b.params <- params

# Extract MPDs over delta transforms, from communities/phylogenies with intraspecific AND seq. error branches
sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
# Separating ultrametric and non-ultrametric results
is.ultrametric(sim.Results[[1]]$phylogenies$orig.phylo)
is.ultrametric(sim.Results[[910]]$phylogenies$orig.phylo)

sim.phylos <- lapply(sim.Results, function(x) x$phylogenies)



lapply(sim.Results[[100]]$phylogenies, is.ultrametric)
lapply(sim.Results[[910]]$phylogenies, is.ultrametric)

null.test <- function(results){
  #browser()
  counter <- 0
  for(i in 1:length(results)){
    if(is.null(results[[i]]$phylogenies$orig.phylo)){
      counter <- (counter + 1)
    }
  }
  return(counter)
}


which(is.ultrametric(sim.Results$phylogenies$orig.phylo))

sim.ultra.data <- Filter(is.ultrametric, sim.data$phylogenies)
sim.nonultra.data <- Filter(Negate(is.ultrametric), sim.data$phylogenies)

length(sim.Results)
length(sim.ultra.data)
length(sim.nonultra.data)

# Extract diversity metrics from untransformed phylogenies, for calculating response metrics
sim.MPDs <- lapply(sim.Results, function(x) x$values$MPDs)

# Build results matrix, from which linear model variables will be pulled
results <- params[rep(1:nrow(params), each=30),]
results$delta <- rep(deltas,length(sim.data))
# Make subset of ultrametric results (death rates = 0)
ultra.results <- results[which(results$comm.death == 0),]
ultra.results <- ultra.results[which(ultra.results$intra.death == 0),]
ultra.results <- ultra.results[which(ultra.results$seq.death == 0),]
# Make subset of non-ultrametric results (death rates != 0)
nonultra.results <- results[which(results$comm.death != 0),]
nonultra.results <- nonultra.results[which(nonultra.results$intra.death != 0),]
nonultra.results <- nonultra.results[which(nonultra.results$seq.death != 0),]
# Create term for diversification (birth rate/death rate) for both intraspecific and seq.err.
nonultra.results$comm.div <- nonultra.results$comm.birth/nonultra.results$comm.death
nonultra.results$intra.div <- nonultra.results$intra.birth/nonultra.results$intra.death
nonultra.results$seq.div <- nonultra.results$seq.birth/nonultra.results$seq.death

sim.Results.DT <- setDT(sim.Results)
sim.Results.DT <- as.data.table(sim.Results)




# %%% Plotting %%%----
# *** MPD correlations ***
# Raw correlation values against delta, intraspecific diversification, and seq. err. diversification
plot(results$correl ~ results$delta)
boxplot(results$correl ~ results$delta)
plot(results$correl ~ results$intra.div, col="red")
boxplot(results$correl ~ results$intra.div)
plot(results$correl ~ results$seq.div, col="blue")
boxplot(results$correl ~ results$seq.div)

# *** Rank shifts ***
# Raw rank shifts values against delta, intraspecific diversification, and seq. err. diversification
plot(results$rank.shifts ~ results$delta)
boxplot(results$rank.shifts ~ results$delta)
plot(results$rank.shifts ~ results$intra.div, col="red")
boxplot(results$rank.shifts ~ results$intra.div)
plot(results$rank.shifts ~ results$seq.div, col="blue")
boxplot(results$rank.shifts ~ results$seq.div)

# *** 6 plots: both response metrics, with median and standard deviation values ***
par(mfcol=c(3,2), mar = c(4,2,0,4), oma = c(0, 4, .5, .5), mgp = c(2, 0.6, 0))

# MPD CORRELATIONS
# Include correl response metric value at delta=1.0 in the results matrix, for each simulation instance. 
comp.correls <- results[which(results$delta == 1),18]
comp.correls <- rep(comp.correls, each=30)
results$comp.correls <- comp.correls

# Calculate median/sd of correl values MINUS values at delta=1.0, for original 
correl.d.medians <- tapply((results$correl-results$comp.correls), results$delta, median)
correl.d.sdevs <- tapply((results$correl-results$comp.correls), results$delta, sd)

plot(correl.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-0.25, 0.25)), ylab="", xlab="Delta")
arrows(x0=unique(results$delta), y0=correl.d.medians-correl.d.sdevs, 
       x1=unique(results$delta), y1=correl.d.medians+correl.d.sdevs, 
       code=3, angle=90, length=0.1)

# Calculate median/sd of correl values MINUS values at delta=1.0, for intraspecific transforms
correl.i.medians <- tapply(results$correl-results$comp.correls, results$intra.div, median)
correl.i.sdevs <- tapply(results$correl-results$comp.correls, results$intra.div, sd)

plot(correl.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-0.25, 0.25)), ylab="", xlab="Intraspecific diversification rate")
arrows(x0=unique(results$intra.div), y0=correl.i.medians-correl.i.sdevs, 
       x1=unique(results$intra.div), y1=correl.i.medians+correl.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

# Calculate median/sd of correl values MINUS values at delta=1.0, for seq. err. transforms
correl.s.medians <- tapply(results$correl-results$comp.correls, results$seq.div, median)
correl.s.sdevs <- tapply(results$correl-results$comp.correls, results$seq.div, sd)

plot(correl.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-0.25, 0.25)), ylab="", xlab="Sequencing error diversification rate")
arrows(x0=unique(results$seq.div), y0=correl.s.medians-correl.s.sdevs, 
       x1=unique(results$seq.div), y1=correl.s.medians+correl.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")
mtext("Median MPD correlations", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.55)

# RANK SHIFTS
# Include response delta=1.0 values in the results matrix, for each simulation instance
comp.rank.shifts <- results[which(results$delta == 1),20]
comp.rank.shifts <- rep(comp.rank.shifts, each=30)
results$comp.rank.shifts <- comp.rank.shifts

# Calculate median/sd of rank.shifts values MINUS values at delta=1.0, for original 
rank.shifts.d.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, median)
rank.shifts.d.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, sd)

plot(rank.shifts.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-5, 5)), ylab="", xlab="Delta")
arrows(x0=unique(results$delta), y0=rank.shifts.d.medians-rank.shifts.d.sdevs, 
       x1=unique(results$delta), y1=rank.shifts.d.medians+rank.shifts.d.sdevs, 
       code=3, angle=90, length=0.1)

# Calculate median/sd of rank.shifts values MINUS values at delta=1.0, for intraspecific transforms 
rank.shifts.i.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$intra.div, median)
rank.shifts.i.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$intra.div, sd)

plot(rank.shifts.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-5, 5)), ylab="", xlab="Intraspecific diversification rate")
arrows(x0=unique(results$intra.div), y0=rank.shifts.i.medians-rank.shifts.i.sdevs, 
       x1=unique(results$intra.div), y1=rank.shifts.i.medians+rank.shifts.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

# Calculate median/sd of rank.shifts values MINUS values at delta=1.0, for seq. err. transforms 
rank.shifts.s.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$seq.div, median)
rank.shifts.s.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$seq.div, sd)

plot(rank.shifts.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-5, 5)), ylab="", xlab="Sequencing error diversification rate")
arrows(x0=unique(results$seq.div), y0=rank.shifts.s.medians-rank.shifts.s.sdevs, 
       x1=unique(results$seq.div), y1=rank.shifts.s.medians+rank.shifts.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")
mtext("Median site ranking shifts", side = 2, outer = TRUE, line = -27, cex = 1.0, adj = 0.55)

# *** Plot raw MPD vs. delta ***
# Capture mean MPD values across sites for each simulation instance
mean.MPDs <- t(sapply(sim.data, function(x) apply(x,2,mean,na.rm="TRUE")))

# Capture mean MPD values across sites for original communities (prior to transformations)
mean.MPDs.Original <- sapply(sim.MPDs, mean)

# Subtract mean.MPDs.Original values from each column of mean.MPDs
final.MPDs <- matrix(nrow=675, ncol=30, dimnames=dimnames(mean.MPDs))
for(i in 1:ncol(mean.MPDs)){
  final.MPDs[,i] <- mean.MPDs[,i] - mean.MPDs.Original
}

# Plot lines representing centered mean MPD values for each simulation iteration versus delta
plot(final.MPDs[1,] ~ unique(results$delta), ylim=c(-20,20), type="n", ylab="Mean MPD Values", xlab="Delta")
for(i in 1:nrow(final.MPDs)){
  lines(unique(results$delta), final.MPDs[i,], col=rgb(red=0.3, green=0.1, blue=0.4, alpha=0.1))
}

# %%% Backup variables %%%----
sim.Results <- backup
params <- b.params
