library(geiger)
library(pez)

# Simulate a tree
tree <- sim.bdtree(100)

# Simulate abundances
species.abundances <- round(t(sim.char(tree, 5, model="BM", nsim=50)[,1,]) * 100)
?sim.char
# "simulating evolution of discrete or continuous characters on a phylogenetic tree". 
# Returns "an array of simulated data, either two or three-dimensional. First dimension is the number of taxa, the second the number of characters, and the third the number of simulated data sets.
# Then, we're transposing the 2nd dimension
# Then, we're multiplying by 100 and rounding
# But...how does evolution of character traits simulate species abundances? Or perhaps this is just a means of obtaining those abundance values?
species.abundances[species.abundances < 0] <- 0
# Remove negative abundances

# Make formatting pretty
colnames(species.abundances) <- tree$tip.label
rownames(species.abundances) <- paste("site_", seq_len(nrow(species.abundances)))
species.abundances

# Make comparative.data object
c.data <- comparative.comm(tree, species.abundances)

# Calculate MPD
.mpd(c.data, abundance.weighted=TRUE)
#...there are others: .mntd, .pd, .psv... I could go on!