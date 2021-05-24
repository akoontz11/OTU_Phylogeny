# ANALYSIS OF SIMULATION OUTPUT 

library(xtable)
library(viridis)

# %%% Worker functions %%%----
# Function to calculate correlations of MPD measurements between each delta value (x)
# and "reference" MPD values (r, taken prior to transformations)
.new.correls <- function(x,r){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],r, use="na.or.complete")
  }
  return(value)
}

# Function to calculate difference in site (MPD) rankings 
.ranking.diff <- function(x,r){
  # Generate a numeric vector of baseline site rankings (i.e. prior to transformations)
  baseline.rank <- rank(r)
  # Generate a vector to capture mean changes between rankings
  meanRankChanges <- vector(length = ncol(x))
  for (i in 1:ncol(x)){
    # Calculate the mean of the absolute changes between two different rankings, and place into vector
    meanRankChanges[i] <- mean(abs(baseline.rank - rank(x[,i])))
  }
  return(meanRankChanges)
}

# %%% Ultrametric runs only %%%----
# %%% Read in simulation data %%%
load("OTU_Phylogeny/simResults/simResults.ULTRA_ONLY.RData")
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

# %%% MPD correlations between site and baseline %%%
# Calculate MPD correlations on simulation data
t.correls <- mapply(.new.correls, sim.data, sim.MPDs)
# Match the correlations to the parameters
results$correl <- as.numeric(t.correls)
# Model effects on MPD correlations
s.model.correl <- lm(correl ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(s.model.correl)

# %%% Ranking differences between site and baseline %%%
# Calculate mean number of site rank shifts
rank.shifts <- mapply(.ranking.diff, sim.data, sim.MPDs)
# Match the ranking differences to the parameters
results$rank.shifts <- as.numeric(rank.shifts)
# Model effects on site diversity rankings
s.model.rankShifts <- lm(rank.shifts ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(s.model.rankShifts)

# %%% Non-ultrametric runs %%%----
# %%% Read in simulation data %%%
load("OTU_Phylogeny/simResults/simResults.RData")
# Generate backup data of simulation variables
backup <- sim.Results
b.params <- params

# NON-ULTRAMETRIC INSTANCES
# Extract MPD values over delta transformations for "final" community/phylogey set
sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
# Subset parameters matrix by successful (i.e. not NULL) simulation instances
params <- params[sapply(sim.data, is.matrix),]
# Remove NAs from simulation data
sim.data <- sim.data[which(!is.na(sim.data))]
# Extract original MPD values of each community, prior to branch additions or transformations
sim.MPDs <- lapply(sim.Results, function(x) x$values$MPDs)
sim.MPDs <- sim.MPDs[sapply(sim.MPDs, Negate(is.null))]
# Build results matrix, from which linear model variables will be pulled
results <- params[rep(1:nrow(params), each=30),]
results$delta <- rep(deltas,length(sim.data))
# Create diversification rate terms
results$intra.div <- results$intra.birth/results$intra.death
results$seq.div <- results$seq.birth/results$seq.death

# %%% MPD correlations between site and baseline %%%
# Calculate MPD correlations on simulation data
t.correls <- mapply(.new.correls, sim.data, sim.MPDs)
results$correl <- as.numeric(t.correls)
# Model effects on MPD correlations
n.model.correl <- lm(correl ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(n.model.correl)

# %%% Ranking differences between site and baseline %%%
# Calculate mean number of site rank shifts
rank.shifts <- mapply(.ranking.diff, sim.data, sim.MPDs)
results$rank.shifts <- as.numeric(rank.shifts)
# Model effects on site diversity rankings
n.model.rankShifts <- lm(rank.shifts ~ scale(log10(delta))+scale(intra.div)+scale(seq.div)+scale(comm.spp),data=results,na.action=na.omit)
summary(n.model.rankShifts)

# ULTRAMETRIC INSTANCES
# Extract MPD values over delta transformations for "final" community/phylogey set
sim.ultra.data <- lapply(sim.ultraResults, function(x) x$transforms$seq.transform)
# Subset parameters matrix by successful (i.e. not NULL) simulation instances
ultra.params <- ultra.params[sapply(sim.ultra.data, is.matrix),]
# Remove NAs from simulation data
sim.ultra.data <- sim.ultra.data[which(!is.na(sim.ultra.data))]
# Extract original MPD values of each community, prior to branch additions or transformations
sim.ultraMPDs <- lapply(sim.ultraResults, function(x) x$values$MPDs)
sim.ultraMPDs <- sim.ultraMPDs[sapply(sim.ultraMPDs, Negate(is.null))]
# Build results matrix, from which linear model variables will be pulled
ultra.results <- ultra.params[rep(1:nrow(ultra.params), each=30),]
ultra.results$delta <- rep(deltas,length(sim.ultra.data))

# %%% MPD correlations between site and baseline %%%
# Calculate MPD correlations on simulation data
u.correls <- mapply(.new.correls, sim.ultra.data, sim.ultraMPDs)
ultra.results$correl <- as.numeric(u.correls)
# Model effects on MPD correlations
u.model.correl <- lm(correl ~ scale(log10(delta))+scale(comm.spp),data=ultra.results,na.action=na.omit)
summary(u.model.correl)

# %%% Ranking differences between site and baseline %%%
# Calculate mean number of site rank shifts
u.rank.shifts <- mapply(.ranking.diff, sim.ultra.data, sim.ultraMPDs)
ultra.results$rank.shifts <- as.numeric(u.rank.shifts)
# Model effects on site diversity rankings
u.model.rankShifts <- lm(rank.shifts ~ scale(log10(delta))+scale(comm.spp),data=ultra.results,na.action=na.omit)
summary(u.model.rankShifts)

# %%% Non-ultrametric, SESmpd %%%----
# %%% Read in simulation data %%%
load("OTU_Phylogeny/simResults/simResults_20210520.RData")
# Generate backup data of simulation variables
nu.backup <- sim.Results 
u.backup <- sim.ultra.Results 
nu.params <- params
u.params <- ultra.params

# NON-ULTRAMETRIC INSTANCES
# Remove params rows of instances in error
params <- params[which(lapply(sim.Results, length) != 1),]
# Remove errored instances
sim.Results <- sim.Results[which(lapply(sim.Results, length) != 1)]

# Extract MPD values over delta transformations for "final" community/phylogey set
sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
# Extract original SESmpd values of each community, prior to branch additions or transformations
sim.SESmpds <- lapply(sim.Results, function(x) x$values$SESmpds)
# Extract SESmpd values of phylogenies with intra+seq branches at delta=1.0
sim.SESmpds.seq <- lapply(sim.Results, function(x) x$transforms$seq.transform[,10])
# Extract phylogenies with intra and seq branches appended 
sim.seq.phylos <- lapply(sim.Results, function(x) x$phylogenies$seq.phylo)

# Build results matrix, from which linear model variables will be pulled
results <- params[rep(1:nrow(params), each=30),]
results$delta <- rep(deltas, length(sim.data))
# Capture (total) number of species 
results$tips <- sapply(sim.seq.phylos, Ntip)
# Create diversification rate terms
# results$inter.div <- results$inter.birth/results$inter.death
results$intra.div <- results$intra.birth/results$intra.death
results$seq.div <- results$seq.birth/results$seq.death

# %%% MPD correlations between site and baseline %%%
# Calculate MPD correlations on simulation data
t.correls <- mapply(.new.correls, sim.data, sim.SESmpds)
results$correl <- as.numeric(t.correls)
# Model effects on MPD correlations
nu.correl.orig <- lm(correl ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                       (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                          scale(intra.div)+scale(seq.div)+scale(tips)),data=results,na.action=na.omit)
summary(nu.correl.orig)
xtable(nu.correl.orig, caption="")

# %%% MPD correlations between site and seq at delta=1.0 %%%
# Calculate MPD correlations on simulation data
seq.correls <- mapply(.new.correls, sim.data, sim.SESmpds.seq)
results$correl.seq <- as.numeric(seq.correls)
# Model effects on MPD correlations
nu.correl.seq <- lm(correl.seq ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                      (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                         scale(intra.div)+scale(seq.div)+scale(tips)),data=results,na.action=na.omit)
summary(nu.correl.seq)
xtable(nu.correl.seq, caption="")

# %%% Ranking differences between site and baseline %%%
# Calculate mean number of site rank shifts
rank.shifts <- mapply(.ranking.diff, sim.data, sim.SESmpds)
results$rank.shifts <- as.numeric(rank.shifts)
# Model effects on site diversity rankings
nu.rankShifts.orig <- lm(rank.shifts ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                           (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                              scale(intra.div)+scale(seq.div)+scale(tips)),data=results,na.action=na.omit)
summary(nu.rankShifts.orig)
xtable(nu.rankShifts.orig, caption="")

# %%% Ranking differences between site and seq at delta=1.0 %%%
# Calculate mean number of site rank shifts
rank.shifts.seq <- mapply(.ranking.diff, sim.data, sim.SESmpds.seq)
results$rank.shifts.seq <- as.numeric(rank.shifts.seq)
# Model effects on site diversity rankings
nu.rankShifts.seq <- lm(rank.shifts.seq ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                          (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                             scale(intra.div)+scale(seq.div)+scale(tips)),data=results,na.action=na.omit)
summary(nu.rankShifts.seq)
xtable(nu.rankShifts.seq, caption="")

# ULTRAMETRIC INSTANCES
# Remove params rows of instances in error
ultra.params <- ultra.params[which(lapply(sim.ultra.Results, length) != 1),]
# Remove errored instances
sim.ultra.Results <- sim.ultra.Results[which(lapply(sim.ultra.Results, length) != 1)]

# Extract MPD values over delta transformations for "final" community/phylogey set
sim.ultra.data <- lapply(sim.ultra.Results, function(x) x$transforms$seq.transform)
# Extract SESmpd values of each community, prior to branch additions or transformations
sim.ultra.SESmpds <- lapply(sim.ultra.Results, function(x) x$values$SESmpds)
# Extract SESmpd values of each community, after branch additions and at delta = 1.0
sim.ultra.SESmpds.seq <- lapply(sim.ultra.Results, function(x) x$transforms$seq.transform[,10])
# Extract phylogenies with intra and seq branches appended 
u.sim.seq.phylos <- lapply(sim.ultra.Results, function(x) x$phylogenies$seq.phylo)

# Build results matrix, from which linear model variables will be pulled
ultra.results <- ultra.params[rep(1:nrow(ultra.params), each=30),]
ultra.results$delta <- rep(deltas, length(sim.ultra.data))
# Capture (total) number of species 
ultra.results$tips <- sapply(u.sim.seq.phylos, Ntip)
# Calculate Pybus and Harvey (2000) gamma statistic
ultra.results$gammas <- sapply(u.sim.seq.phylos, gammaStat)

# %%% MPD correlations between site and baseline %%%
# Calculate MPD correlations on simulation data
u.correls <- mapply(.new.correls, sim.ultra.data, sim.ultra.SESmpds)
ultra.results$correl <- as.numeric(u.correls)
# Model effects on MPD correlations (all terms: initial species, birth rates, tips, gamma)
u.correl.orig <- lm(correl ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                      (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                         scale(intra.birth)+scale(seq.birth)+scale(tips)+scale(gammas)),
                    data=ultra.results,na.action=na.omit)
summary(u.correl.orig)
xtable(u.correl.orig, caption="")

# %%% MPD correlations between site and seq at delta=1.0 %%%
# Calculate MPD correlations on simulation data
u.correls.seq <- mapply(.new.correls, sim.ultra.data, sim.ultra.SESmpds.seq)
ultra.results$correl.seq <- as.numeric(u.correls.seq)
# Model effects on MPD correlations (all terms: initial species, birth rates, tips, gamma)
u.correl.seq <- lm(correl.seq ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                     (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                        scale(intra.birth)+scale(seq.birth)+scale(tips)+scale(gammas)),
                   data=ultra.results,na.action=na.omit)
summary(u.correl.seq)
xtable(u.correl.seq, caption="")

# %%% Ranking differences between site and baseline %%%
# Calculate mean number of site rank shifts
u.rank.shifts <- mapply(.ranking.diff, sim.ultra.data, sim.ultra.SESmpds)
ultra.results$rank.shifts <- as.numeric(u.rank.shifts)
# Model effects on site diversity rankings (all terms: initial species, birth rates, tips, gamma)
u.rankShifts.orig <- lm(rank.shifts ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                          (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                             scale(intra.birth)+scale(seq.birth)+scale(tips)+scale(gammas)),
                        data=ultra.results,na.action=na.omit)
summary(u.rankShifts.orig)
xtable(u.rankShifts.orig, caption="")

# %%% Ranking differences between site and seq at delta=1.0 %%%
# Calculate mean number of site rank shifts
u.rank.shifts.seq <- mapply(.ranking.diff, sim.ultra.data, sim.ultra.SESmpds.seq)
ultra.results$rank.shifts.seq <- as.numeric(u.rank.shifts.seq)
# Model effects on site diversity rankings (all terms: initial species, birth rates, tips, gamma)
u.rankShifts.seq <- lm(rank.shifts.seq ~ (I(scale(log10(delta))^2)+scale(log10(delta)))*
                         (scale(intra.ratio)+scale(seq.ratio)+scale(comm.spp)+
                            scale(intra.birth)+scale(seq.birth)+scale(tips)+scale(gammas)),
                       data=ultra.results,na.action=na.omit)
summary(u.rankShifts.seq)
xtable(u.rankShifts.seq, caption="")

# %%% Plotting--Old %%%----
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
# %%% SESmpd Analysis Plots %%%----
# 1 row, 2 columns
par(mfcol=c(2,1))

# SESmpd values of seq trees over delta
# Plot SESmpd means for each site over delta values--non-ultrametric
# ymin <- min(sapply(sim.data, min)); ymax <- max(sapply(sim.data, max))
for(i in 1:length(sim.data)){
  means <- apply(sim.data[[i]], 2, mean)
  if(i == 1){
    plot(means ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric", ylim=c(-5,5), pch=20, col=i)
    lines(deltas, means, col=i)
  } else {
    points(means ~ deltas, col=i, pch=20)
    lines(deltas, means, col=i)
  }
}

# Plot SESmpd means for each site over delta values--ultrametric
# ymin <- min(sapply(sim.ultra.data, min)); ymax <- max(sapply(sim.ultra.data, max))
for(i in 1:length(sim.ultra.data)){
  means <- apply(sim.ultra.data[[i]], 2, mean)
  if(i == 1){
    plot(means ~ deltas, xlab="Delta", ylab="SESmpd", main="Ultrametric", ylim=c(-5,5), pch=20, col=i)
    lines(deltas, means, col=i)
  } else {
    points(means ~ deltas, col=i, pch=20)
    lines(deltas, means, col=i)
  }
}

# Capture mean SESmpd values across sites for non-ultrametric instances
nu.mean.SESmpds <- t(sapply(sim.data, function(x) apply(x,2,mean,na.rm="TRUE")))
# Plot raw values using transparencies (alpha)
plot(nu.mean.SESmpds[1,] ~ unique(results$delta), ylim=c(-4,4), 
     type="n", ylab="Mean SESmpd Values", xlab="Delta", main="Non-ultrametric SESmpd")
for(i in 2:nrow(nu.mean.SESmpds)){
  lines(unique(results$delta), nu.mean.SESmpds[i,], 
        col=rgb(red=0.3, green=0.1, blue=0.4, alpha=0.1))
}

# Capture mean SESmpd values across sites for ultrametric instances
u.mean.SESmpds <- t(sapply(sim.ultra.data, function(x) apply(x,2,mean,na.rm="TRUE")))
# Plot raw values using transparencies (alpha)
plot(u.mean.SESmpds[1,] ~ unique(ultra.results$delta), ylim=c(-4,4), 
     type="n", ylab="Mean SESmpd Values", xlab="Delta", main="Ultrametric SESmpd")
for(i in 2:nrow(u.mean.SESmpds)){
  lines(unique(ultra.results$delta), u.mean.SESmpds[i,], 
        col=rgb(red=0.3, green=0.1, blue=0.4, alpha=0.1))
}

# # Troubleshooting
# plot(apply(sim.data[[25]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 25", ylim=c(-1.5,1.5), pch=20, col=4)
# lines(deltas, apply(sim.data[[25]], 2, mean), col=4)
# phy.d.transform.plot(sim.Results[[25]]$phylogenies$seq.phylo, sim.Results[[25]]$abundances$seq.community, deltas, "Sim Instance #25")
# 
# plot(apply(sim.data[[15]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 15", ylim=c(-3.5,1.5), pch=20, col=4)
# lines(deltas, apply(sim.data[[15]], 2, mean), col=4)
# phy.d.transform.plot(sim.Results[[15]]$phylogenies$seq.phylo, sim.Results[[15]]$abundances$seq.community, deltas, "Sim Instance #15")
# 
# which(sapply(sim.seq.phylos, Ntip) < 300)
# sim.issues <- sim.Results[which(sapply(sim.seq.phylos, Ntip) < 300)]
# sim.issue.data <- lapply(sim.issues, function(x) x$transforms$seq.transform)
# 
# for(i in 1:length(sim.issue.data)){
#   means <- apply(sim.issue.data[[i]], 2, mean)
#   if(i == 1){
#     plot(means ~ deltas, xlab="Delta", ylab="SESmpd", main="", ylim=c(-1.5,1.5), pch=20, col=i)
#     lines(deltas, means, col=i)
#   } else {
#     points(means ~ deltas, col=i, pch=20)
#     lines(deltas, means, col=i)
#   }
# }
# 
# plot(apply(sim.data[[1]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 1", ylim=c(-1.5,1.5), pch=20, col=3)
# lines(deltas, apply(sim.data[[1]], 2, mean), col=3)
# 
# plot(apply(sim.data[[73]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 1", ylim=c(-1.5,1.5), pch=20, col=3)
# lines(deltas, apply(sim.data[[73]], 2, mean), col=3)
# 
# phy.d.transform.plot(sim.Results[[15]]$phylogenies$seq.phylo, sim.Results[[15]]$abundances$seq.community, deltas, "Sites correspond with one another")
# phy.d.transform.plot(sim.Results[[73]]$phylogenies$seq.phylo, sim.Results[[73]]$abundances$seq.community, deltas, "Trends differ between sites")
# 
# # No change in SESmpd: site 16, 25
# plot(apply(sim.data[[16]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 16", ylim=c(ymin,ymax), pch=20, col=3)
# lines(deltas, apply(sim.data[[16]], 2, mean), col=3)
# phy.d.transform.plot(sim.Results[[16]]$phylogenies$seq.phylo, sim.Results[[16]]$abundances$seq.community, deltas, "Sim Instance #16")
# params[16,]
# 
# plot(sim.Results[[16]]$phylogenies$orig.phylo, show.tip.label = F)
# plot(sim.Results[[16]]$phylogenies$intra.phylo, show.tip.label = F)
# plot(sim.Results[[16]]$phylogenies$seq.phylo, show.tip.label = F)
# Ntip(sim.Results[[16]]$phylogenies$seq.phylo)
# 
# plot(apply(sim.data[[25]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 25", ylim=c(ymin,ymax), pch=20, col=3)
# lines(deltas, apply(sim.data[[25]], 2, mean), col=3)
# phy.d.transform.plot(sim.Results[[25]]$phylogenies$seq.phylo, sim.Results[[25]]$abundances$seq.community, deltas, "Sim Instance #25")
# params[25,]
# 
# plot(sim.Results[[25]]$phylogenies$orig.phylo, show.tip.label = F)
# plot(sim.Results[[25]]$phylogenies$intra.phylo, show.tip.label = F)
# plot(sim.Results[[25]]$phylogenies$seq.phylo, show.tip.label = F)
# Ntip(sim.Results[[25]]$phylogenies$seq.phylo)
# 
# # Appreciable change in SESmpd: site 15, 32
# plot(apply(sim.data[[15]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 15", ylim=c(ymin,ymax), pch=20, col=6)
# lines(deltas, apply(sim.data[[15]], 2, mean), col=6)
# phy.d.transform.plot(sim.Results[[15]]$phylogenies$seq.phylo, sim.Results[[15]]$abundances$seq.community, deltas, "Sim Instance #15")
# params[15,]
# 
# plot(sim.Results[[15]]$phylogenies$orig.phylo, show.tip.label = F)
# plot(sim.Results[[15]]$phylogenies$intra.phylo, show.tip.label = F)
# plot(sim.Results[[15]]$phylogenies$seq.phylo, show.tip.label = F)
# Ntip(sim.Results[[15]]$phylogenies$seq.phylo)
# 
# plot(apply(sim.data[[32]], 2, mean) ~ deltas, xlab="Delta", ylab="SESmpd", main="Non-ultrametric, instance 32", ylim=c(ymin,ymax), pch=20, col=6)
# lines(deltas, apply(sim.data[[32]], 2, mean), col=6)
# phy.d.transform.plot(sim.Results[[32]]$phylogenies$seq.phylo, sim.Results[[32]]$abundances$seq.community, deltas, "Sim Instance #32")
# params[32,]
# 
# plot(sim.Results[[32]]$phylogenies$orig.phylo, show.tip.label = F)
# plot(sim.Results[[32]]$phylogenies$intra.phylo, show.tip.label = F)
# plot(sim.Results[[32]]$phylogenies$seq.phylo, show.tip.label = F)
# Ntip(sim.Results[[32]]$phylogenies$seq.phylo)

# 2 rows, 2 columns
par(mfcol=c(2,2), mar = c(2,2,2,2), oma = c(0, 4, .5, .5), mgp = c(2, 0.6, 0))

# Correlations between original SESmpd metrics and baseline
# Non-ultrametric
correls <- mapply(.new.correls, sim.data, sim.SESmpds)
ymin <- (min(correls)-1.5); ymax <- (max(correls)+1.5)
for(i in 1:length(sim.data)){
  if(i == 1){
    plot(correls[,i] ~ deltas, xlab="Delta", ylab="SESmpd correlations", 
         main="", ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, correls[,i], col=i)
  } else {
    points(correls[,i] ~ deltas, col=i, pch=20)
    lines(deltas, correls[,i], col=i)
  }
}

# Ultrametric
u.correls <- mapply(.new.correls, sim.ultra.data, sim.ultra.SESmpds)
ymin <- (min(u.correls)-1.5); ymax <- (max(u.correls)+1.5)
for(i in 1:length(sim.ultra.data)){
  if(i == 1){
    plot(u.correls[,i] ~ deltas, xlab="Delta", ylab="SESmpd correlations", 
         main="", ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, u.correls[,i], col=i)
  } else {
    points(u.correls[,i] ~ deltas, col=i, pch=20)
    lines(deltas, u.correls[,i], col=i)
  }
}

# Site shifts in SESmpd ranking between sites at all delta values and baseline
# Non-ultrametric
rank.shifts <- mapply(.ranking.diff, sim.data, sim.SESmpds)
ymin <- (min(rank.shifts)-1.5); ymax <- (max(rank.shifts)+1.5)
for(i in 1:length(sim.data)){
  if(i == 1){
    plot(rank.shifts[,i] ~ deltas, main="", xlab="", ylab="",
         ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, rank.shifts[,i], col=i)
  } else {
    points(rank.shifts[,i] ~ deltas, col=i, pch=20)
    lines(deltas, rank.shifts[,i], col=i)
  }
}

# Ultrametric
u.rank.shifts <- mapply(.ranking.diff, sim.ultra.data, sim.ultra.SESmpds)
ymin <- (min(u.rank.shifts)-1.5); ymax <- (max(u.rank.shifts)+1.5)
for(i in 1:length(sim.ultra.data)){
  if(i == 1){
    plot(u.rank.shifts[,i] ~ deltas, main="", xlab="", ylab="",
         ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, u.rank.shifts[,i], col=i)
  } else {
    points(u.rank.shifts[,i] ~ deltas, col=i, pch=20)
    lines(deltas, u.rank.shifts[,i], col=i)
  }
}

mtext("Non-ultrametric", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.85)
mtext("Ultrametric", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.25)
mtext("Correlations", side = 3, outer = TRUE, line = -1, cex = 1.0, adj = 0.20)
mtext("Rank shifts", side = 3, outer = TRUE, line = -1, cex = 1.0, adj = 0.80)

# 2 rows, 2 columns
par(mfcol=c(2,2), mar = c(2,2,2,2), oma = c(0, 4, .5, .5), mgp = c(2, 0.6, 0))

# Correlations between original SESmpd metrics and delta=1.0
# Non-ultrametric
correls <- mapply(.new.correls, sim.data, sim.SESmpds.seq)
ymin <- (min(correls)-1.5); ymax <- (max(correls)+1.5)
for(i in 1:length(sim.data)){
  if(i == 1){
    plot(correls[,i] ~ deltas, xlab="Delta", ylab="SESmpd correlations", 
         main="", ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, correls[,i], col=i)
  } else {
    points(correls[,i] ~ deltas, col=i, pch=20)
    lines(deltas, correls[,i], col=i)
  }
}

# Ultrametric
u.correls <- mapply(.new.correls, sim.ultra.data, sim.ultra.SESmpds.seq)
ymin <- (min(u.correls)-1.5); ymax <- (max(u.correls)+1.5)
for(i in 1:length(sim.ultra.data)){
  if(i == 1){
    plot(u.correls[,i] ~ deltas, xlab="Delta", ylab="SESmpd correlations", 
         main="", ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, u.correls[,i], col=i)
  } else {
    points(u.correls[,i] ~ deltas, col=i, pch=20)
    lines(deltas, u.correls[,i], col=i)
  }
}

# Site shifts in SESmpd ranking between sites at all delta values and delta=1.0
# Non-ultrametric
rank.shifts <- mapply(.ranking.diff, sim.data, sim.SESmpds.seq)
ymin <- (min(rank.shifts)-1.5); ymax <- (max(rank.shifts)+1.5)
for(i in 1:length(sim.data)){
 if(i == 1){
   plot(rank.shifts[,i] ~ deltas, main="", xlab="", ylab="",
        ylim=c(ymin,ymax), pch=20, col=i)
   lines(deltas, rank.shifts[,i], col=i)
 } else {
   points(rank.shifts[,i] ~ deltas, col=i, pch=20)
   lines(deltas, rank.shifts[,i], col=i)
 }
}

# Ultrametric
u.rank.shifts <- mapply(.ranking.diff, sim.ultra.data, sim.ultra.SESmpds.seq)
ymin <- (min(u.rank.shifts)-1.5); ymax <- (max(u.rank.shifts)+1.5)
for(i in 1:length(sim.ultra.data)){
  if(i == 1){
    plot(u.rank.shifts[,i] ~ deltas, main="", xlab="", ylab="",
         ylim=c(ymin,ymax), pch=20, col=i)
    lines(deltas, u.rank.shifts[,i], col=i)
  } else {
    points(u.rank.shifts[,i] ~ deltas, col=i, pch=20)
    lines(deltas, u.rank.shifts[,i], col=i)
  }
}

mtext("Non-ultrametric", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.85)
mtext("Ultrametric", side = 2, outer = TRUE, line = 2, cex = 1.0, adj = 0.25)
mtext("Correlations", side = 3, outer = TRUE, line = -1, cex = 1.0, adj = 0.20)
mtext("Rank shifts", side = 3, outer = TRUE, line = -1, cex = 1.0, adj = 0.80)
