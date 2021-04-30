# ADD BRANCHES TO PHYLOGENETIC TREE (INTRASPECIFIC/SEQUENCING ERROR)

library(ape)
library(geiger)
library(phytools)

# New function, rescaling tree----
# Append trees onto phylogeny tips. birth, death, steps are fed to sim.bdtree. type determines populations or seq. error
add.branch <- function(tree,birth,death,type=c("pops","seq.err")){
  type <- match.arg(type)
  # If argument `tree` is without branches (i.e. NULL), return a tree that is NULL. Otherwise, add branches to it.
  if(is.null(tree)){
    return(tree)
  }else{
    tips <- Ntip(tree)
    # Make copy of original tree, to keep track of the tip names on the original tree
    # Also, get longest branch of original tree (for rescaling later)
    orig.tree <- tree
    old.max.depth <- max(branching.times(orig.tree))
    for(i in 1:tips){
      orig.tip <- orig.tree$tip.label[i]
      # Generate a new tree to be added onto the tip of the original tree, based on arguments
      new.tree <- sim.bdtree(b=birth, d=death, stop="time", t=1)
      # Rename new tree tip labels
      # Get number of tips of new tree
      new.tips <- Ntip(new.tree)
      # For each tip in new tree, get the current tip, and reformat its name
      for (j in 1:new.tips){
        cur.tip <- new.tree$tip.label[j]
        # The `type` argument denotes whether to add a "p" or an "e" to the new tips, for populations or seq. error.
        if(type=="seq.err"){
          cur.tip <- sub("s","e",cur.tip)
        } else{
          cur.tip <- sub("s","p",cur.tip)
        }
        new.name <- paste0(orig.tip,"--",cur.tip)
        # Add that new name to the tip of the new tree
        new.tree$tip.label[j] <- new.name
      }
      tree <- bind.tree(tree,new.tree,where=1)
      # (where always equals 1, because with every addition, next "open" tip becomes tip number 1)
    }
    # Get longest branch of new tree
    new.max.depth <- max(branching.times(tree))
    # To maintain interspecific to intra/seq ratios, new tree must be one branch length unit longer than old tree
    if(new.max.depth != (old.max.depth+1)){
      # If not, scale the new tree to make it so
      tree$edge.length <- (tree$edge.length)*((old.max.depth+1)/new.max.depth)
    }
    return(tree)
  }
}

# # Old function----
# add.branch <- function(tree,birth,death,steps,type=c("pops","seq.err")){
#   type <- match.arg(type)
#   # If argument `tree` is without branches (i.e. NULL), return a tree that is NULL. Otherwise, add branches to it.
#   if(is.null(tree)){
#     return(tree)
#   }else{
#     tips <- Ntip(tree)
#     # Make copy of original tree, to keep track of the tip names on the original tree.
#     orig.tree <- tree
#     for(i in 1:tips){
#       orig.tip <- orig.tree$tip.label[i]
#       # Generate a new tree to be added onto the tip of the original tree, based on arguments.
#       # new.tree <- sim.bdtree(b=birth, d=death, stop="time", t=steps)
#       new.tree <- sim.bdtree(b=birth, d=death, stop="time", t=steps)
#       # Rename new tree tip labels.
#       # Get number of tips of new tree.
#       new.tips <- Ntip(new.tree)
#       # For each tip in new tree, get the current tip, and reformat its name.
#       for (j in 1:new.tips){
#         cur.tip <- new.tree$tip.label[j]
#         # The `type` argument denotes whether to add a "p" or an "e" to the new tips, for populations or seq. error.
#         if(type=="seq.err"){
#           cur.tip <- sub("s","e",cur.tip)
#         } else{
#           cur.tip <- sub("s","p",cur.tip)
#         }
#         new.name <- paste0(orig.tip,"--",cur.tip)
#         # Add that new name to the tip of the new tree.
#         new.tree$tip.label[j] <- new.name
#       }
#       tree <- bind.tree(tree,new.tree,where=1)
#       # (where always equals 1, because with every addition, next "open" tip becomes tip number 1.)
#     }
#     # Drop extinct branches of the tree, such that only extant leaves are shown
#     # Compare the number of extinct tips to the total number of tips
#     # if(Ntip(tree)-length(is.extinct(tree))<=1){
#     #   # If tree has one or less branches remaining, return NULL for the resulting tree
#     #   tree <- NULL
#     #   }else{
#     #     # Otherwise, prune the tree of tips that have gone extinct
#     #     tree <- drop.extinct(tree)
#     #     }
#     return(tree)
#   }
# }