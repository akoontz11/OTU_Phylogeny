# Script for reading in Bonnie's microbial dataset

library(ape)
setwd("/home/akoontz11/OTU_Phylogeny/Microbial_Dataset/")

# Read in and reformat community matrices
# Bacterial
bact.comm <- as.matrix(read.csv2("CRbact_communitydat.csv", header=TRUE, sep=","))
soil.cores <- bact.comm[,1]
bact.comm <- apply(bact.comm[,-1],2,as.numeric)
rownames(bact.comm) <- soil.cores
sort(unique(c(bact.comm)))
# Fungal
fung.comm <- as.matrix(read.csv2("CRfungi_communitydat.csv", header=TRUE, sep=","))
soil.cores <- fung.comm[,1]
fung.comm <- apply(fung.comm[,-1],2,as.numeric)
rownames(fung.comm) <- soil.cores
colnames(fung.comm) <- gsub(pattern="X", replacement="", colnames(fung.comm))
sort(unique(c(fung.comm)))

# Read in phylogenetic trees (aligned with mafft, built with RAxML)
# Bacterial
# ...
# Fungal
fung.phylo <- read.tree("phylos/fungal/RAxML_bestTree.repCRfung_Feb28_phy")
# Prune tree of tips not present in the community matrix
fung.phylo <- keep.tip(fung.phylo, colnames(fung.comm))
plot(fung.phylo, show.tip.label=T)




