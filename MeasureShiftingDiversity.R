# %%% ANALYSIS OF SIMULATION OUTPUT %%%

setwd("/home/akoontz11/OTU_Phylogeny/")
library(xtable)
library(viridis)

# %%% FUNCTION DECLARATIONS %%% ----
# Function for standardizing variables in linear models
z.transform <- function(data){
  s.data <- (data-mean(data))/sd(data)
  return(s.data)
}

# A worker function that will calculate the correlations between each delta value (x)
# and a "reference" delta value (d, taken prior to transformations)
.new.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d, use="na.or.complete")
  }
  return(value)
}

# A worker function that will calculate the difference in site rankings 
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

# Worker function for comparing two mpd values: prior to transformation, and at delta=1.0 after intra and seq additions
.vertical.comparison <- function(x,d){
  # x is vector of mpd values measured after delta transformations; d is original mpd value prior to transformations
  # The use argument allows correlations to be caluculated off of all pairs that are present
  value <- cor(x[,10], d, use="na.or.complete")
  return(value)
}

# %%% READING IN AND FORMATTING OF SIMULATION DATA %%% ----
# Pared down simulation results (51 MB; excluding instances of death rates > birth rates)
load("simResults.RData")
# Generate backup data of simulation variables
backup <- sim.Results
b.params <- params

# Removed erred lines from output (to be addressed later...)
params <- params[sapply(sim.Results, Negate(is.character)),]
sim.Results <- Filter(Negate(is.character), sim.Results)

# Extract relevant data from the results
#sim.data <- lapply(sim.Results, function(x) x$transforms$orig.transform)
#sim.data <- lapply(sim.Results, function(x) x$transforms$intra.transform)
sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)

# Remove non-matrix simulation values
params <- params[sapply(sim.data, is.matrix),]
sim.data <- Filter(is.matrix, sim.data)

# Extract diversity metrics from untransformed phylogenies (for later comparison)
sim.MPDs <- lapply(sim.Results, function(x) x$values$MPDs)

# Build results matrix, from which linear model variables will be pulled
results <- params[rep(1:nrow(params), each=30),]
results$delta <- rep(deltas,length(sim.data))
# Create term for diversification rate (birth rate/death rate) for both intraspecific and sequening error differences
results$intra.div <- results$intra.birth/results$intra.death
results$seq.div <- results$seq.birth/results$seq.death

# %%% CORRELATIONS OF MPD VALUES BETWEEN SITE AND BASELINE %%% ----
# Using mapply on worker function calculating correlations between delta values
t.correls <- mapply(.new.correls, sim.data, sim.MPDs)
# Match the correlations to the parameters
results$correl <- as.numeric(t.correls)
# Include delta=1.0 values in the results matrix, for each simulation instance
comp.correls <- results[which(results$delta == 1),18]
comp.correls <- rep(comp.correls, each=30)
results$comp.correls <- comp.correls

# %%% LINEAR MODELS AND SUMMARIES %%%
# Model terms: log10(delta), intra diversification, seq diversification, number of initial species
s.model.correl <- lm(correl ~ z.transform(log10(delta))+z.transform(intra.div)+z.transform(seq.div)+z.transform(comm.spp),data=results,na.action=na.omit)
summary(s.model.correl)
xtable(s.model.correl)

# %%% PLOTTING %%%
# Plotting raw correl values
# delta
plot(results$correl ~ results$delta)
boxplot(results$correl ~ results$delta)
# intra.div
plot(results$correl ~ results$intra.div, col="red")
boxplot(results$correl ~ results$intra.div)
# seq.div
plot(results$correl ~ results$seq.div, col="blue")
boxplot(results$correl ~ results$seq.div)

# Plotting medians and standard deviations
# Plotting the medians of correl values minus correl values at delta=1.0
correl.d.medians <- tapply((results$correl-results$comp.correls), results$delta, median)
# Standard deviations of correl values minus correl values at delta = 1.0
correl.d.sdevs <- tapply((results$correl-results$comp.correls), results$delta, sd)

plot(correl.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-0.25, 0.25)))
arrows(x0=unique(results$delta), y0=correl.d.medians-correl.d.sdevs, 
       x1=unique(results$delta), y1=correl.d.medians+correl.d.sdevs, 
       code=3, angle=90, length=0.1)

# Plotting the centered medians of correl values grouped by intraspecific diversification
correl.i.medians <- tapply(results$correl-results$comp.correls, results$intra.div, median)
# Standard deviations of centered correl values grouped by intraspecific diversification
correl.i.sdevs <- tapply(results$correl-results$comp.correls, results$intra.div, sd)

plot(correl.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-0.25, 0.25)))
arrows(x0=unique(results$intra.div), y0=correl.i.medians-correl.i.sdevs, 
       x1=unique(results$intra.div), y1=correl.i.medians+correl.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

# Plotting the medians of correl values grouped by seq. error diversification
correl.s.medians <- tapply(results$correl-results$comp.correls, results$seq.div, median)
# Standard deviations of correl values grouped by seq. error diversification
correl.s.sdevs <- tapply(results$correl-results$comp.correls, results$seq.div, sd)

plot(correl.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-0.25, 0.25)))
arrows(x0=unique(results$seq.div), y0=correl.s.medians-correl.s.sdevs, 
       x1=unique(results$seq.div), y1=correl.s.medians+correl.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")

# %%% DIFFERENCE IN SITE RANKINGS (I.E. CROSSINGS OVER) BETWEEN SITE AND BASELINE %%% ----
# Using mapply on worker function determining number of site rank shiftings
rank.shifts <- mapply(.ranking.diff, sim.data, sim.MPDs)
# Match the ranking differences to the parameters
results$rank.shifts <- as.numeric(rank.shifts)
# Include delta=1.0 values in the results matrix, for each simulation instance
comp.rank.shifts <- results[which(results$delta == 1),20]
comp.rank.shifts <- rep(comp.rank.shifts, each=30)
results$comp.rank.shifts <- comp.rank.shifts

# %%% LINEAR MODELS AND SUMMARIES %%%
# Model terms: log10(delta), intra diversification, seq diversification, number of initial species
s.model.rankShifts <- lm(rank.shifts ~ z.transform(log10(delta))+z.transform(intra.div)+z.transform(seq.div)+z.transform(comm.spp),data=results,na.action=na.omit)
summary(s.model.rankShifts)
xtable(s.model.rankShifts)

# %%% PLOTTING %%%
# Plotting raw rank.shifts values
# delta
plot(results$rank.shifts ~ results$delta)
boxplot(results$rank.shifts ~ results$delta)
# intra.div
plot(results$rank.shifts ~ results$intra.div, col="red")
boxplot(results$rank.shifts ~ results$intra.div)
# seq.div
plot(results$rank.shifts ~ results$seq.div, col="blue")
boxplot(results$rank.shifts ~ results$seq.div)

# Plotting medians and standard deviations
# Plotting the medians of rank.shifts values minus rank.shifts values at delta=1.0
rank.shifts.d.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, median)
# Standard deviations of rank.shifts values minus rank.shifts values at delta = 1.0
rank.shifts.d.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, sd)

plot(rank.shifts.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-5, 5)))
arrows(x0=unique(results$delta), y0=rank.shifts.d.medians-rank.shifts.d.sdevs, 
       x1=unique(results$delta), y1=rank.shifts.d.medians+rank.shifts.d.sdevs, 
       code=3, angle=90, length=0.1)

# Plotting the centered medians of rank.shifts values grouped by intraspecific diversification
rank.shifts.i.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$intra.div, median)
# Standard deviations of centered rank.shifts values grouped by intraspecific diversification
rank.shifts.i.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$intra.div, sd)

plot(rank.shifts.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-5, 5)))
arrows(x0=unique(results$intra.div), y0=rank.shifts.i.medians-rank.shifts.i.sdevs, 
       x1=unique(results$intra.div), y1=rank.shifts.i.medians+rank.shifts.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

# Plotting the centered medians of rank.shifts values grouped by seq. error diversification
rank.shifts.s.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$seq.div, median)
# Standard deviations of centered rank.shifts values grouped by seq. error diversification
rank.shifts.s.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$seq.div, sd)

plot(rank.shifts.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-5, 5)))
arrows(x0=unique(results$seq.div), y0=rank.shifts.s.medians-rank.shifts.s.sdevs, 
       x1=unique(results$seq.div), y1=rank.shifts.s.medians+rank.shifts.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")

# Four plots: both response variables, grouped both by delta and by simulated error
par(mfcol=c(3,2), mar = c(4,2,0,4), oma = c(0, 4, .5, .5), mgp = c(2, 0.6, 0))

plot(correl.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-0.25, 0.25)), ylab="", xlab="Delta")
arrows(x0=unique(results$delta), y0=correl.d.medians-correl.d.sdevs, 
       x1=unique(results$delta), y1=correl.d.medians+correl.d.sdevs, 
       code=3, angle=90, length=0.1)

plot(correl.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-0.25, 0.25)), ylab="", xlab="Intraspecific diversification rate")
arrows(x0=unique(results$intra.div), y0=correl.i.medians-correl.i.sdevs, 
       x1=unique(results$intra.div), y1=correl.i.medians+correl.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

plot(correl.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-0.25, 0.25)), ylab="", xlab="Sequencing error diversification rate")
arrows(x0=unique(results$seq.div), y0=correl.s.medians-correl.s.sdevs, 
       x1=unique(results$seq.div), y1=correl.s.medians+correl.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")

mtext("Median MPD correlations", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.55)

plot(rank.shifts.d.medians ~ unique(results$delta), pch=16, ylim=range(c(-5, 5)), ylab="", xlab="Delta")
arrows(x0=unique(results$delta), y0=rank.shifts.d.medians-rank.shifts.d.sdevs, 
       x1=unique(results$delta), y1=rank.shifts.d.medians+rank.shifts.d.sdevs, 
       code=3, angle=90, length=0.1)

plot(rank.shifts.i.medians ~ unique(results$intra.div), pch=16, col="red", ylim=range(c(-5, 5)), ylab="", xlab="Intraspecific diversification rate")
arrows(x0=unique(results$intra.div), y0=rank.shifts.i.medians-rank.shifts.i.sdevs, 
       x1=unique(results$intra.div), y1=rank.shifts.i.medians+rank.shifts.i.sdevs, 
       code=3, angle=90, length=0.1, col="red")

plot(rank.shifts.s.medians ~ unique(results$seq.div), pch=16, col="blue", ylim=range(c(-5, 5)), ylab="", xlab="Sequencing error diversification rate")
arrows(x0=unique(results$seq.div), y0=rank.shifts.s.medians-rank.shifts.s.sdevs, 
       x1=unique(results$seq.div), y1=rank.shifts.s.medians+rank.shifts.s.sdevs, 
       code=3, angle=90, length=0.1, col="blue")
mtext("Median site ranking shifts", side = 2, outer = TRUE, line = -27, cex = 1.0, adj = 0.55)

# %%% PLOTTING RAW MPD VERSUS DELTA %%% ----
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

# %%% COMPARISON OF ORIGINAL MPD VALUES TO VALUES AFTER BRANCH ADDITION %%% ----
# Using mapply on worker function determining number of site rank shiftings
v.comps <- mapply(.vertical.comparison, sim.data, sim.MPDs)
# Match the ranking differences to the parameters
results$v.comps <- as.numeric(v.comps)

# %%% LINEAR MODELS AND SUMMARIES %%%
model.verticalComparison <- lm(v.comps ~ z.transform(intra.div)+z.transform(seq.div), data=results, na.action=na.omit)
summary(model.verticalComparison)

# %%% BACKUP VARIABLES %%% ----
sim.Results <- backup
params <- b.params
