# Adding branches to phylogenetic tree (intraspecific/sequencing error)

# Load packages
library(ape)
library(geiger)
library(phytools)

# Function for appending trees onto the tips of an existing tree. 
# The birth, death, and steps arguments all describe the parameters of the simulated trees to be appended to the original tree.
add.branch <- function(tree,birth,death,steps,type=c("pops","seq.err")){
  # Set type of branch addition according to argument (this specifies only the letters used for the tips of the appended tree)
  type <- match.arg(type)
  # Determine number of species present in original tree. Additional trees will be added to each originally present tip.
  tips <- Ntip(tree)
  # Make copy of original tree. This is in order to keep track of the tip names on the original tree.
  orig.tree <- tree
  # For each species in original tree,
  for(i in 1:tips){
    # get name of original tip.
    orig.tip <- orig.tree$tip.label[i]
    # Generate a new tree to be added onto the tip of the original tree, based on the function arguments specified.
    new.tree <- sim.bdtree(b=birth, d=death, stop="time", t=steps)
    # Rename tip labels of new tree to be prefixed by tip label on the original tree, where new tree will be bound.
    # Get number of tips of new tree.
    new.tips <- Ntip(new.tree)
    # For each tip in new tree,
    for (j in 1:new.tips){
     # get the current tip, and reformat its name.
     cur.tip <- new.tree$tip.label[j]
     # The seq.err flag denotes whether to add a "p" or an "e" to the new tips
     # (This is meant to represent either populations or sequencing error being added to the tree, respectively)
     if(type=="seq.err"){
       cur.tip <- sub("s","e",cur.tip) 
     } else{
       cur.tip <- sub("s","p",cur.tip)
     }
     new.name <- paste0(orig.tip,"--",cur.tip)
     # Add that new name to the tip of the new tree.
     new.tree$tip.label[j] <- new.name
    }
    # Bind newly generated tree to current tip of original tree.
    tree <- bind.tree(tree,new.tree,where=1)
    # The where argument of this command always equals 1, because every time a new tree is appended to the original, 
    # the next "open" tip of the original tree (which is where the next new tree will be added) becomes tip number 1.
  }
  # Drop extinct branches of the tree, such that only extant leaves are shown
  tree <- drop.extinct(tree)
  return(tree)
}

tree <- sim.bdtree(b=0.1, d=0, stop="time", t=10)
plot(tree, main="Original tree")

intra.tree <- add.branch(tree,birth=0.1,death=0.2,steps=5,type="pops")
plot(intra.tree, main="Population tree")

seq.tree <- add.branch(intra.tree,birth=0.2,death=0.2,steps=3,type="seq.err")
plot(seq.tree, main="Sequencing error tree")
