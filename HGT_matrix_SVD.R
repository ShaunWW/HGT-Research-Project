#This Script Makes a HGT matrix for the SVD split data

setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits/svd_with_linebreaks")
# ^ SVD little splits


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
# > combo.vec[1,]
# [1] 1 2 3 4


support.matrix <- matrix(data=0, nrow= 1001, ncol = 4)

x<-12
while( x < 1002)
{
  # Read Tree File
  quart <- combo.vec[x,]
  if ( !is.element(13, quart))
  {
    split.file <- paste(c("littlesplits_",quart[1],"_",quart[2],"_", quart[3],"_",quart[4],".txt"), sep='', collapse = "")
    split.table <-read.table(split.file, sep = "", stringsAsFactors = FALSE)
    L <- length(split.table[[1]])
    AB.sup <- 0
    AC.sup <- 0
    AD.sup <- 0
    # For each orthologog set in a given quartet, make a vector of those which
    # show more than 80% support
    for(ortholog in 1:L)
    {
      if (split.table[ortholog,2] > 79)
      {
        # We want to save and count all splits supporting a given quartet
        split.num <- split.table[ortholog,1]
        # split2 will hold the info to update matrix
        if (split.num== 1 ){ AB.sup <- AB.sup +1 }
        else if ( split.num == 2 ){ AC.sup <- AC.sup +1}
        else if ( split.num == 3 ){ AD.sup <- AD.sup +1}
        
      }
    }
    
    # update support matrix:
    support.matrix[x,1]<- paste(quart, sep = "_", collapse = "_")
    support.matrix[x,2]<- AB.sup
    support.matrix[x,3]<- AC.sup
    support.matrix[x,4]<- AD.sup
  }
  x <- x+1
}

write.table(support.matrix, file = "SVD_Support_data.txt")
#^^^^^
# We now have Support Matrix, Check which subsplits show support for HGT

Likelihood<- function(n1, n2, n3, p1, p2, p3){
  
  P = (p1^n1)*(p2^n2)*(p3^n3)
  #(factorial(n1+n2+n3)/(factorial(n1)*factorial(n2)*factorial(n3)))*  gets cancelled by division
  return(P)
}
support.matrix <- read.table(file = "SVd_Support_data", stringsAsFactors = FALSE)
fam.table<- support.matrix

SVD.Mat <- matrix( data = 0, nrow = 14, ncol = 14)

for (quart in 1:1001)
{
  # Check for Null entries
  if (fam.table[quart,2]+fam.table[quart,3]+fam.table[quart,4] !=0)
  {
    vals <- sort(c(fam.table[quart,2],fam.table[quart,3],fam.table[quart,4]))
    mean <- (vals[1] + vals[2])/2
    tot <- sum(vals)
    #Mean could be zero (100 0 0 split)
    if (mean != 0)
    {
      # Run Likelihood test
      Ptest <- 0.5 *((vals[1]+vals[2])/100) # corsponds to p of vals1&2, vals3 is the max
      top <- Likelihood(vals[1], vals[2], vals[3], Ptest , Ptest, vals[3]/tot)
      bot <- Likelihood(vals[1], vals[2], vals[3], vals[1]/tot , vals[2]/tot, vals[3]/tot)
      
      Pval <- -2*log(top/bot)
      
      # Check Pval vs Threshold
      if(Pval > 3.841)
      {
        # If this loop is entered:
        # Ratio test implys points are different
        # Need to find "second" split
        x <- 2
        while(x < 5)
        {
          if (fam.table[quart,x] == vals[2])
          {
            split.num <- x
          }
          x<- x +1
        }
        
        # split list for tracking splits
        # split.num.list[quart]<- split.num
        
        species <- strsplit(fam.table[quart,1], split = "_")
        # split2 will hold the info to update matrix
        if (split.num == 2 )
        {
          split2a <-c(species[[1]][1], species[[1]][2])
          split2b <-c(species[[1]][3], species[[1]][4])
        }else if ( split.num ==3 )
        {
          split2a <-c(species[[1]][1], species[[1]][3])
          split2b <-c(species[[1]][2], species[[1]][4])
        }else if ( split.num==4 )
        {
          split2a <-c(species[[1]][1], species[[1]][4])
          split2b <-c(species[[1]][2], species[[1]][3])
        }
        
        # Record Splits that Show Support
        rowa<- as.numeric(split2a[1])
        cola<- as.numeric(split2a[2])
        rowb<- as.numeric(split2b[1])
        colb<- as.numeric(split2b[2])
        
        SVD.Mat[rowa, cola] <- SVD.Mat[rowa, cola] +1
        SVD.Mat[rowb, colb] <- SVD.Mat[rowb, colb] +1
      } 
    }
  }
}

write.table(SVD.Mat, file = "SVD_Mat_95.txt", sep = "    ")