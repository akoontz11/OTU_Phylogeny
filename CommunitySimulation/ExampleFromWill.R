# Make the parent matrix
mat <- matrix(1:25, nrow=5)
colnames(mat) <- letters[1:5]

# Assume we know how many each row will expand into; pre-allocate matrix
expansions <- 1:5
big.mat <- matrix(NA, ncol=sum(expansions), nrow=nrow(mat))

# Do work looping over sites
for(i in seq_len(nrow(mat))){
  # Setup tracker index (x)
  x <- 1
  # Loop over the true species (columns)
  for(j in seq_len(ncol(mat))){
    # Make a random shuffle of the appropriate size (note the unname it matters!)
    rnd.shuffle <- sample(x=c(unname(mat[i,j]), rep(0,expansions[j]-1)), size=expansions[j])
    # Put new values in
    cat("Inner loop iteration: ",j,"\n")  
    print(big.mat)
    big.mat[i,seq(x,length.out=expansions[j])] <- rnd.shuffle
    # Increment tracker
    x <- x+expansions[j]
  }
}
mat
big.mat



expansions
expansions[4]

i=1
j=2
rnd.shuffle <- sample(x=c(unname(mat[i,j]), rep(0,expansions[j]-1)), size=expansions[j])
# x is a concatenation of the matrix value and some series of zeroes (length equal to expansion size-1)
# We take 'size' samples from x, equal to expansions value. This makes up the rnd.shuffle vector. 
rnd.shuffle

x=4
j=3
seq(x,length.out=expansions[j])

i=1
j=1
x=1

unname(mat[1,4])
rep(0,expansions[j]-1)

rnd.shuffle <- sample(x=c(unname(mat[i,j]), rep(0,expansions[j]-1)), size=expansions[j])
rnd.shuffle

c(unname(mat[i,j]), rep(0,expansions[j]-1))

big.mat[i,seq(x,length.out=expansions[j])] <- rnd.shuffle
big.mat
rnd.shuffle

seq_len(nrow(mat))
?seq_len
?unname
?sample
