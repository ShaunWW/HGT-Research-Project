# This script reads the most supported quartete split for the 1001 combos of tree_data
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits")

# Generate combo.vec
L<- 14
combo.vec<- matrix(data=0, nrow= 1001, ncol = 4)
a<-1
b<-2
c<-3
d<-4
index<-1
while (a <= (L-3))
{
  while( b <= (L-2))
  {
    while(c <= (L-1))
    {
      while(d <= L)
      {
        combo.vec[index ,] <- c(a,b,c,d)
        
        index<- index+1
        d<- d+1
      }
      c<- c+1
      d<- c+1
    }
    c<- b+1
    b<- b+1
  }
  b<- a+1
  a<- a+1
}
rm(a,b,c,d, index, L)

count.vector <- matrix( data = 0, nrow = 1001, ncol = 4)
# For every TreeData.file
for ( P in 1:1970)
{
  # Read in the tree.data
  split.table <- read.table(file = paste(c("TreeData_", P ,".txt"), sep ="", collapse = ""), stringsAsFactors = FALSE)
  support.table <- matrix( data = 0, nrow = 1001, ncol = 3)
  
  #For each embedded Quaret
  for ( i in 1:1001)
  {
    # read quart Name
    quart<-c(combo.vec[i,1],combo.vec[i,2], combo.vec[i,3], combo.vec[i,4])
    support.table[i,1]<- paste(quart, sep = "", collapse = "_")
    
    # Find entry >80
    s<-2
    support <- FALSE
    while( s < 5 && support == FALSE)
    {
      if (split.table[i,s] >= 80)
      {
        support <- TRUE
        entry <- s
      }
      s <- s+1
    }
    
    # Record Entry or no strong support
    if( support == TRUE)
    { # ^ Support is high for a given split
     
      # Record split numbers & update count vector
      if (entry ==2 )
      {# AB split
        support.table[i,2] <- paste(c(quart[1], quart[2]), sep ="", collapse = "_")
        count.vector[i,2] <- as.numeric(count.vector[i,2]) +1
      } else if (entry == 3)
      {
        support.table[i,2] <- paste(c(quart[1], quart[3]), sep ="", collapse = "_")
        count.vector[i,3] <- as.numeric(count.vector[i,3]) +1
      } else
      {
        support.table[i,2] <- paste(c(quart[1], quart[4]), sep ="", collapse = "_")
        count.vector[i,4] <- as.numeric(count.vector[i,4]) +1
      }
      
      # Record support value
      support.table[i,3]<- split.table[i, entry]
    } else # No strong support for quart
    {
      support.table[i,2]<- "No Strong Support"
      # support.table[i,3] <- "N/A
      
      # Record Highest level of support
      support.table[i,3]<- max( split.table[i,2], split.table[i,3], split.table[i,4])
    }
  }
  # Write tree.data to file
  # write.table(support.table, file = paste(c("SplitSupport_", P ,".txt"), sep ="", collapse = ""), sep = "   ")
  P <- P +1  
}

# set quartet names
count.vector[,1]<- support.table[,1]

# Write Support Table
write.table(count.vector, file = "Metafile_90.txt" , sep = "   ")


#count.vector <- read.table(file = "Metafile.txt", stringsAsFactors = FALSE)
quart.support <- matrix(data = 0, nrow = 1001, ncol = 2)
quart.support[,1]<- count.vector[,1]
# find most supported split for each quartet
for (j in 1:1001)
{
  AB <- count.vector[j,2]
  AC <- count.vector[j,3]
  AD <- count.vector[j,4]
  
  cherry <- "0"
  high.sup <- max(AB,AC,AD)
  match<-0 # use match to check for ties
  if ( high.sup == AB)
  {
    cherry <- "AB"
    match <- match +1
  }
  if ( high.sup == AC)
  {
    cherry <- "AC"
    match <- match +1
  }
  if ( high.sup == AD)
  {
    cherry <- "AD"
    match <- match +1
  }
  
  if( match < 2)
  { #No tie was found
    quart.support[j,2]<- cherry 
  }
  
}

write.table(quart.support , file = "Quartet Best Split_90.txt", sep = "    ")
