# %%% Analysis of simulation output %%%

# Function for standardizing
z.transform <- function(data){
  s.data <- (data-mean(data))/sd(data)
  return(s.data)
}

# Load simulation results, and make a copy of the relevant output variable
load("demoresults.RData")
# Generate backup data of simulation variables
backup <- sim.Results
b.params <- params
#--------------------------------------------------------------------------
# Removed erred lines from output (to be addressed later...)
params <- params[sapply(sim.Results, Negate(is.character)),]
sim.Results <- Filter(Negate(is.character), sim.Results)

# Extract relevant data from the results
#data <- lapply(sim.Results, function(x) x$transforms$orig.transform)
data <- lapply(sim.Results, function(x) x$transforms$intra.transform)
#data <- lapply(sim.Results, function(x) x$transforms$seq.transform)
params <- params[sapply(data, is.matrix),]
data <- Filter(is.matrix, data)

# ---CORRELATIONS OF MPD VALUES BETWEEN SITE AND BASELINE---
# A worker function that will calculate the correlations between each delta value and "baseline" delta value (1.0)
.site.correls <- function(x){
  value <- numeric(length=ncol(x))
  for (i in 1:ncol(x)){
    value[i] <- cor(x[,i],x[,10])
  }
  return(value)
}
# Calculate the correlations
correls <- sapply(data, .site.correls)
# Match the correlations to the parameters
results <- params[rep(1:nrow(params), each=30),]
results$correl <- as.numeric(correls)
results$delta <- rep(deltas,length(data))

# Generating model using standardized variables
s.model.correl <- lm(correl ~ z.transform(delta)+I(z.transform(delta)^2)+z.transform(intra.birth)+z.transform(intra.death)+z.transform(intra.steps)+z.transform(seq.birth)+z.transform(seq.death)+z.transform(seq.steps),data=results)
summary(s.model.correl)
# Plotting command
#plot(results$correl~results$delta)

# ---DIFFERENCE IN SITE RANKINGS (I.E. CROSSINGS OVER) BETWEEN SITE AND BASELINE---
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
rank.shifts <- sapply(data, .ranking.diff)
# Match the ranking differences to the parameters
results <- params[rep(1:nrow(params), each=30),]
results$rank.shifts <- as.numeric(rank.shifts)
results$delta <- rep(deltas,length(data))

# Generating model using standardized variables
s.model.rank.shifts <- lm(rank.shifts ~ z.transform(delta)+I(z.transform(delta)^2)+z.transform(intra.birth)+z.transform(intra.death)+z.transform(intra.steps)+z.transform(seq.birth)+z.transform(seq.death)+z.transform(seq.steps),data=results)
summary(s.model.rank.shifts)
# Plotting command
#plot(results$rank.shifts~results$delta)
#plot(results$rank.shifts~results$intra.birth)

sim.Results <- backup
params <- b.params

