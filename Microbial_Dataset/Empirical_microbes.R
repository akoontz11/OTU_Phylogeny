# Script for reading in Bonnie's microbial dataset

library(ape)
#setwd("/home/akoontz11/OTU_Phylogeny/Microbial_Dataset/")

bact.comm <- as.matrix(read.csv2("CRbact_communitydat.csv", header=TRUE, sep=","))
soil.cores <- bact.comm[,1]
bact.comm <- apply(bact.comm[,-1],2,as.numeric)
rownames(bact.comm) <- soil.cores
sort(unique(c(bact.comm)))

bact.phylo <- read.tree("CRbact_repset.tre")
plot(bact.phylo)

