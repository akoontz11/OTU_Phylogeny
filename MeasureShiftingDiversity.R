# %%% ANALYSIS OF SIMULATION OUTPUT %%%

setwd("/home/akoontz11/OTU_Phylogeny/")
library(xtable)

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
# # Orig.transforms: delta, (delta)^2
# o.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2), data=results, na.action=na.omit)
# summary(o.model.correl)
# 
# # Intra.transforms: delta, (delta)^2, intra diversification
# i.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div), data=results, na.action=na.omit)
# summary(i.model.correl)

# All model terms
s.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div)+z.transform(seq.div),data=results,na.action=na.omit)
summary(s.model.correl)
xtable(s.model.correl)

# %%% PLOTTING %%%
# par(mar=c(5.0, 4.0, 1.0, 2.0) + 0.1)
# 
# plot(tapply(s.model.correl$model$correl, results$delta, mean) ~ unique(results$delta), 
#      ylab="Model terms: correlations", xlab="Delta Value", pch=16)
# 
# plot(tapply(s.model.correl$model$correl, results$intra.div, mean) ~ unique(results$intra.div),
#      ylab="Model terms: correlations", xlab="Intraspecific Diversification", pch=16, col="red")
# 
# plot(tapply(s.model.correl$model$correl, results$seq.div, mean) ~ unique(results$seq.div), 
#      ylab="Model terms: correlations", xlab="Sequencing Error Diversification", pch=16, col="blue")
# 
# # Three plots in one window
# par(mfcol=c(3,1), mar=c(3, 4, 1, 8) + 0.1)
# 
# plot(tapply(s.model.correl$model$correl, results$delta, mean) ~ unique(results$delta), 
#      ylab="Model terms: correlations", xlab="Delta Value", pch=16)
# 
# plot(tapply(s.model.correl$model$correl, results$intra.div, mean) ~ unique(results$intra.div),
#       ylab="Model terms: correlations", xlab="Intraspecific Diversification", pch=16, col="red")
# 
# plot(tapply(s.model.correl$model$correl, results$seq.div, mean) ~ unique(results$seq.div), 
#      ylab="Model terms: correlations", xlab="Sequencing Error Diversification", pch=16, col="blue")

# UPDATED PLOTS
# Plotting raw data
plot(results$correl ~ results$delta)
boxplot(results$correl ~ results$delta)

# Plotting the medians of correl values minus correl values at delta=1.0, with error bars
correl.medians <- tapply((results$correl-results$comp.correls), results$delta, median)
# Standard deviations of correl values minus correl values at delta = 1.0
correl.c.sdevs <- tapply((results$correl-results$comp.correls), results$delta, sd)
# Standard deviations of raw correl values
#correl.sdevs <-  apply(t.correls, 1, sd)

plot(correl.medians ~ unique(results$delta), pch=16, ylim=range(c(-1.0, 1.0)))

arrows(x0=unique(results$delta), y0=correl.medians-correl.c.sdevs, 
       x1=unique(results$delta), y1=correl.medians+correl.c.sdevs, 
       code=3, angle=90, length=0.1)

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
# # Orig.transforms: delta, (delta)^2
# o.model.rankShifts <- lm(rank.shifts ~ z.transform(log(delta))+I(z.transform(log(delta))^2), data=results, na.action=na.omit)
# summary(o.model.rankShifts)
# 
# # Intra.transforms: delta, (delta)^2, intra diversification
# i.model.rankShifts <- lm(rank.shifts ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div), data=results, na.action=na.omit)
# summary(i.model.rankShifts)

# Seq.transforms: delta, (delta)^2, intra diversification, seq diversification
s.model.rankShifts <- lm(rank.shifts ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div)+z.transform(seq.div),data=results,na.action=na.omit)
summary(s.model.rankShifts)
xtable(s.model.rankShifts)

# %%% PLOTTING %%%
# plot(tapply(s.model.rankShifts$model$rank.shifts, results$delta, mean) ~ unique(results$delta), 
#      ylab="Model terms: ranking shifts", xlab="Delta", pch=16)
# 
# plot(tapply(s.model.rankShifts$model$rank.shifts, results$intra.div, mean) ~ unique(results$intra.div), 
#      ylab="Model terms: ranking shifts", xlab="Intraspecific Diversification", pch=16, col="red")
# 
# plot(tapply(s.model.rankShifts$model$rank.shifts, results$seq.div, mean) ~ unique(results$seq.div), 
#      ylab="Model terms: ranking shifts", xlab="Sequencing Error Diversification", pch=16, col="blue")

# UPDATED PLOTS
# Plotting raw data
plot(results$rank.shifts ~ results$delta)

boxplot(results$rank.shifts ~ results$delta)

# Plotting the medians of rank.shifts values minus rank.shifts values at delta=1.0, with error bars
rank.shifts.medians <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, median)
# Standard deviations of rank.shifts values minus rank.shifts values at delta = 1.0
rank.shifts.c.sdevs <- tapply((results$rank.shifts-results$comp.rank.shifts), results$delta, sd)
# Standard deviations of raw correl values
#rank.shifts.sdevs <-  apply(rank.shifts, 1, sd)

plot(rank.shifts.medians ~ unique(results$delta), pch=16, ylim=range(c(-10, 10)))

arrows(x0=unique(results$delta), y0=rank.shifts.medians-rank.shifts.c.sdevs, 
       x1=unique(results$delta), y1=rank.shifts.medians+rank.shifts.c.sdevs, 
       code=3, angle=90, length=0.1)

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
