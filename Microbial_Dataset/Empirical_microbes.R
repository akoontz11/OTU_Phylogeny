# Script for reading in Bonnie's microbial dataset

library(ape)
#setwd("/home/akoontz11/OTU_Phylogeny/Microbial_Dataset/")

# Read in community matrices
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
sort(unique(c(fung.comm)))

# Read in phylogenetic trees (aligned with mafft, built with RAxML)
# Bacterial
# ...
# Fungal
fung.phylo <- read.tree("phylos/fungal/RAxML_bestTree.repCRfung_Feb28_phy")
plot(fung.phylo, show.tip.label=F)
 
# %%% FUNCTION FOR TESTING TRANSFORMATIONS OVER MICROBIAL COMMUNITIES %%%
MicrobialCommnunity <- function(comm, phy, intra.birth,intra.death,intra.steps,seq.birth,seq.death,seq.steps){
  # %%% CLEANUP COMMUNITY, AND SETUP COMMUNITY WITH POPULATIONS %%%
  # Drop sites with one or less species in them from community matrix  (by removing rows with only one nonzero value in them)
  # (We do this because our mpd calculation is abundance weighted, and mpd can only be calculated between different species)
  comm <- comm[apply(comm, 1, function(x) sum(x!=0) > 1),,drop=F]
  # Drop extinct species from community matrix (by removing columns with only 0s in them)
  comm <- comm[,apply(comm,2,function(x) !all(x==0))]
  # Get the names of remaining species in community (after removing extinctions and 'singletons')
  species.names <- colnames(comm)
   # Trim phylogeny to only contain remaining species (using drop.tip)
   # (Tips are determined by looking for the tip.labels that are not in the species names)
   phy <- drop.tip(phy, tip = phy$tip.label[!phy$tip.label %in% species.names])
   # Get the site names
   sites <- rownames(comm)
   
  # INTRASPECIFIC
  # Add intraspecific difference to simulated community phylogeny
  intra.phy <- add.branch(phy,birth=intra.birth,death=intra.death,steps=intra.steps,"pops")
  # Create an empty matrix for the intra.phy
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  intra.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(intra.phy),dimnames = list(sites,intra.phy$tip.label)),error=function(cond) {return(NULL)})
  population.names <- colnames(intra.comm)
  # Map abundance values from original community matrix to columns of intra.comm 
  intra.comm <- abundance.mapping(species.names, population.names, sites, comm, intra.comm)
  
  # SEQ. ERROR
  # Add random error to intra.tree phylogeny
  seq.phy <- add.branch(intra.phy,birth=seq.birth,death=seq.death,steps=seq.steps,"seq.err")
  # Create an empty matrix for the seq.phy
  # (tryCatch statement included to account for trees in which all species have gone extinct)
  seq.comm <- tryCatch(matrix(nrow=length(sites),ncol=Ntip(seq.phy),dimnames = list(sites,seq.phy$tip.label)),error=function(cond) {return(NULL)})
  individual.names <- colnames(seq.comm)
  # Map abundance values from original community matrix to columns of seq.comm
  seq.comm <- abundance.mapping(species.names, individual.names, sites, comm, seq.comm)
  
  # CAPTURE MPD VALUES OVER TRANSFORMATIONS FOR EACH COMMUNITY TYPE
  # Build mpd versus delta transformation values matrix for original simulated community
  orig.test <- phy.d.transform(phy, comm, deltas)
  # Build mpd versus delta transformation values matrix for community with appended populations
  intra.test <- phy.d.transform(intra.phy, intra.comm, deltas)
  # Build mpd versus delta transformation values matrix for community with sequencing error added on
  seq.test <- phy.d.transform(seq.phy, seq.comm, deltas)
  
  # CAPTURE PHYLOGENETIC DIVERSITY METRICS OF ORIGINAL PHYLOGENIES (FOR LATER COMPARISON)
  # Mean pairwise distance
  orig.MPD <- .mpd(comparative.comm(phy, comm, force.root = 0), abundance.weighted=TRUE)
  # Changing names to match format from phy.d.transform function
  names(orig.MPD) <- paste("Site",1:nrow(comm),sep="_")
  
  # PACKAGING SIMULATION DATA
  # Export a list containing all data, for each community "type" (original, intra, and seq)
  # Abundances
  community.abundances <- list(orig.community=comm,intra.community=intra.comm,seq.community=seq.comm)
  # Phylogenies
  community.phylogenies <- list(orig.phylo=phy,intra.phylo=intra.phy,seq.phylo=seq.phy)
  # MPD matrices, from delta transforms
  community.transforms <- list(orig.transform=orig.test,intra.transform=intra.test,seq.transform=seq.test)
  # Phylogenetic diversity metrics, of original (untransformed) communities/phylogenies
  original.diversityMetrics <- list(MPDs=orig.MPD)
  # Return data
  simulation.data <- list(phylogenies=community.phylogenies,abundances=community.abundances,transforms=community.transforms,values=original.diversityMetrics)
  return(simulation.data)
}