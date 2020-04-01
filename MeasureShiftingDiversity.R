# %%% ANALYSIS OF SIMULATION OUTPUT %%%

setwd("/home/akoontz11/OTU_Phylogeny/")

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
    value[i] <- cor(x[,i],d)
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

# %%% LINEAR MODELS AND SUMMARIES %%%
# Orig.transforms: delta, (delta)^2
# o.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2), data=results, na.action=na.omit)
# summary(o.model.correl)
# 
# # Intra.transforms: delta, (delta)^2, intra diversification
# i.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div), data=results, na.action=na.omit)
# summary(i.model.correl)

# Seq.transforms: delta, (delta)^2, intra diversification, seq diversification
s.model.correl <- lm(correl ~ z.transform(log(delta))+I(z.transform(log(delta))^2)+z.transform(intra.div)+z.transform(seq.div),data=results,na.action=na.omit)
summary(s.model.correl)

# %%% DIFFERENCE IN SITE RANKINGS (I.E. CROSSINGS OVER) BETWEEN SITE AND BASELINE %%% ----
# Using mapply on worker function determining number of site rank shiftings
rank.shifts <- mapply(.ranking.diff, sim.data, sim.MPDs)
# Match the ranking differences to the parameters
results$rank.shifts <- as.numeric(rank.shifts)

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

# %%% OLD WORKER FUNCTIONS %%% ----

# A worker function that will calculate the correlations between each delta value and "baseline" delta value (1.0)
# .site.correls <- function(x){
#   value <- numeric(length=ncol(x))
#   for (i in 1:ncol(x)){
#     value[i] <- cor(x[,i],x[,10])
#   }
#   return(value)
# }
# # Calculate the correlations
# correls <- sapply(sim.data, .site.correls)

# # A worker function that will calculate the difference in site rankings
# .ranking.diff <- function(x){
#   # Generate a numeric vector of baseline site rankings (i.e. delta=1.0)
#   baseline.rank <- as.numeric(rank(sort(x[,10],decreasing = F)))
#   # Generate a vector of names of these site rankings
#   baseline.names <- names(sort(x[,10],decreasing = F))
#   # Generate a vector to capture mean changes between rankings
#   meanRankChanges <- vector(length = ncol(x))
#   for (i in 1:ncol(x)){
#     # Get the names of sites sorted from lowest to highest for the current column (delta value)
#     comparison <- names(sort(x[,i],decreasing = F))
#     # Generate a numeric vector corresponding to how site rankings have changed from "baseline" delta value
#     shifts <- match(baseline.names, comparison)
#     # Calculate absolute changes between two different rankings
#     changes <- abs(shifts - baseline.rank)
#     # Calculate the mean of the absolute changes, and place into vector
#     meanRankChanges[i] <- mean(changes)
#   }
#   return(meanRankChanges)
# }
# # Calculate the ranking differences
# rank.shifts <- sapply(sim.data, .ranking.diff)

