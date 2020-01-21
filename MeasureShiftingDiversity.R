# %%% Analysis of simulation output %%%

# Function for standardizing
z.transform <- function(data){
  s.data <- (data-mean(data))/sd(data)
  return(s.data)
}

# Load simulation results
#load("demoresults.RData")
# Load less bulky, pared down simulation results 
load("simResults.RData")
# Generate backup data of simulation variables
backup <- sim.Results
b.params <- params
# Removed erred lines from output (to be addressed later...)
params <- params[sapply(sim.Results, Negate(is.character)),]
sim.Results <- Filter(Negate(is.character), sim.Results)

# Extract relevant data from the results
sim.data <- lapply(sim.Results, function(x) x$transforms$orig.transform)
#sim.data <- lapply(sim.Results, function(x) x$transforms$intra.transform)
#sim.data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
params <- params[sapply(sim.data, is.matrix),]
sim.data <- Filter(is.matrix, sim.data)
# ----CORRELATIONS OF MPD VALUES BETWEEN SITE AND BASELINE----
# A worker function that will calculate the correlations between each delta value and "baseline" delta value (1.0)
.site.correls <- function(x){
  browser()
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],x[,10])
  }
  return(value)
}
# Development code----
# A worker function that will calculate the correlations between each delta value (x)
# and a "reference" delta value (d, taken prior to transformations)
.new.correls <- function(x,d){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],d)
  }
  return(value)
}
untransformedMPDs <- sim.Results[[1]]$values$MPDs

# Testing new worker function
test.correls <- .new.correls(sim.data[[1]],untransformedMPDs)
test.correls
# Using sapply on new worker function
t.correls <- sapply(sim.data, .new.correls, d=sim.Results$values$MPDs)
# Throws error: "error in cor(x[, i], d) : incompatible dimensions"

# Calculate the correlations----
correls <- sapply(sim.data, .site.correls)
# Match the correlations to the parameters
results <- params[rep(1:nrow(params), each=30),]
results$correl <- as.numeric(correls)
results$delta <- rep(deltas,length(sim.data))

# Generating model using standardized variables
s.model.correl <- lm(correl ~ z.transform(delta)+I(z.transform(delta)^2)+z.transform(intra.birth)+z.transform(intra.death)+z.transform(seq.birth)+z.transform(seq.death),data=results,na.action=na.omit)
#s.model.correl <- lm(correl ~ z.transform(delta)+I(z.transform(delta)^2)+z.transform(intra.birth)+z.transform(intra.death)+z.transform(intra.steps)+z.transform(seq.birth)+z.transform(seq.death)+z.transform(seq.steps),data=results,na.action=na.omit)
summary(s.model.correl)
# Plotting command
#plot(results$correl~results$delta)

# ----DIFFERENCE IN SITE RANKINGS (I.E. CROSSINGS OVER) BETWEEN SITE AND BASELINE----
# A worker function that will calculate the difference in site rankings 
.ranking.diff <- function(x){
  # Generate a numeric vector of baseline site rankings (i.e. delta=1.0)
  baseline.rank <- as.numeric(rank(sort(x[,10],decreasing = F)))
  # Generate a vector of names of these site rankings 
  baseline.names <- names(sort(x[,10],decreasing = F))
  # Generate a vector to capture mean changes between rankings
  meanRankChanges <- vector(length = ncol(x))
  for (i in 1:ncol(x)){
    # Get the names of sites sorted from lowest to highest for the current column (delta value)
    comparison <- names(sort(x[,i],decreasing = F))
    # Generate a numeric vector corresponding to how site rankings have changed from "baseline" delta value
    shifts <- match(baseline.names, comparison)
    # Calculate absolute changes between two different rankings
    changes <- abs(shifts - baseline.rank)
    # Calculate the mean of the absolute changes, and place into vector
    meanRankChanges[i] <- mean(changes)
  }
  return(meanRankChanges)
}
# Calculate the ranking differences
rank.shifts <- sapply(sim.data, .ranking.diff)
# Match the ranking differences to the parameters
results <- params[rep(1:nrow(params), each=30),]
results$rank.shifts <- as.numeric(rank.shifts)
results$delta <- rep(deltas,length(sim.data))

# Generating model using standardized variables
s.model.rank.shifts <- lm(rank.shifts ~ z.transform(delta)+I(z.transform(delta)^2)+z.transform(intra.birth)+z.transform(intra.death)+z.transform(intra.steps)+z.transform(seq.birth)+z.transform(seq.death)+z.transform(seq.steps),data=results)
summary(s.model.rank.shifts)
# Plotting command
#plot(results$rank.shifts~results$delta)
#plot(results$rank.shifts~results$intra.birth)

sim.Results <- backup
params <- b.params

