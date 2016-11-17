# This program will Generate a HGT candidate matrix 
# using both coeff of variation Cv, and variance Sig2

setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits")
# ^ treeData

# Maximum likelihood function:
  # Maximul Likelihood 
Likelihood<- function(n1, n2, n3, p1, p2, p3){
  
  P = (p1^n1)*(p2^n2)*(p3^n3)
  #(factorial(n1+n2+n3)/(factorial(n1)*factorial(n2)*factorial(n3)))*  gets cancelled by division
  return(P)
}

# EXAMPLE: Likelihood(20,30,50,0.25, 0.25, 0.5)/Likelihood(20,30,50,0.2, 0.3, 0.5)


# Generate Empty Matricies 
Ratio.Mat <- matrix( data = 0, nrow = 14, ncol = 14)
Chi.Mat <- Ratio.Mat

#For ortholog family table
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits")
fam.table <- read.table(file = "Metafile_80.txt", stringsAsFactors = FALSE)
  for( quart in 1:1001)
  {
    
    # Check if elements are 0
    if (fam.table[quart,2]+fam.table[quart,3]+fam.table[quart,4] !=0)
    {
      ## If not 0, find two smallest values
      vals <- sort(c(fam.table[quart,2],fam.table[quart,3],fam.table[quart,4]))
        # ^ first two entries are the smaller values
    
      ## Calculate Cv and Sig2 for these values
      mean <- (vals[1] + vals[2])/2
      # Mean could be zero ( 100 0 0 split)
      if (mean != 0)
      {
        ####### Likilhood Test:
        Ptest <- 0.5 *((vals[1]+vals[2])/100) # corsponds to p of vals1&2, vals3 is the max
        top <- Likelihood(vals[1], vals[2], vals[3], Ptest , Ptest, vals[3]/100)
        bot <- Likelihood(vals[1], vals[2], vals[3], vals[1]/100 , vals[2]/100, vals[3]/100)
        
        Pval <- -2*log(top/bot)
        
        ####### GoF
        Chi2  <- ((mean - vals[1])^2)/mean +((mean - vals[2])^2)/mean
        

        # Check Tests against thresholds
        # chi2 for d=1 threshold is 2.7 (10% confidence), Pval follows Chi distribution
        if(Chi2 > 2.7)
        { 
          # If this loop is entered:
          # Chi test implys points are different
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
        
          #split list for tracking splits
          #split.num.list[quart]<- split.num
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
         
            rowa<- as.numeric(split2a[1])
            cola<- as.numeric(split2a[2])
            rowb<- as.numeric(split2b[1])
            colb<- as.numeric(split2b[2])
       
            Chi.Mat[rowa, cola] <- Chi.Mat[rowa, cola] +1
            Chi.Mat[rowb, colb] <- Chi.Mat[rowb, colb] +1
         }

        if(Pval > 2.7)
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

          rowa<- as.numeric(split2a[1])
          cola<- as.numeric(split2a[2])
          rowb<- as.numeric(split2b[1])
          colb<- as.numeric(split2b[2])

          Ratio.Mat[rowa, cola] <- Ratio.Mat[rowa, cola] +1
          Ratio.Mat[rowb, colb] <- Ratio.Mat[rowb, colb] +1
        }
      }
    }
  }
test <- Ratio.Mat-Chi.Mat


      


#Dummy<- Chi.Mat
#Chi.Mat<- Dummy
for (x in 1:14)
{
  # Opposite elements are the same
  for (z in 1:14)
  {
    Chi.Mat[z,x]<- Chi.Mat[x,z]
    Ratio.Mat[z,x] <- Ratio.Mat[x,z]
  }
}

Chi.Tot <- 0
Ratio.Tot <- 0
for (x in 1:14)
{
  # Opposite elements are the same
  for (z in 1:14)
  {
    Chi.Tot<- Chi.Tot + Chi.Mat[x,z]
    Ratio.Tot <- Ratio.Tot + Ratio.Mat[x,z]
  }
}
Chi.Ave <- Chi.Tot/(14*13)
Ratio.Ave <- Ratio.Tot/(14*13)


Chi.Support <- round( Chi.Mat/Chi.Ave, 3)
Ratio.Support <- round( Ratio.Mat/Ratio.Ave , 3)

Agreement.Support <- (Chi.Support+Ratio.Support)/2

for (x in 1:14)
{
  # Diagonal elements are Null
  Chi.Support[x,x] <- "--"
  Ratio.Support[x,x] <- "--"
}



######__________________________________________________________#####
# Write Tables
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits/HGT Support")
# Write Chi.Mat
write.table(Chi.Mat, file = "ChiSquared HGT support.txt", sep = "     ")
write.table(Ratio.Mat, file = "Max Likelihood HGT support.txt", sep = "     ")


# Write Agreement Table:
write.table(Chi.Support, file = "ChiSquared HGT Average.txt", sep = "     ")
write.table(Ratio.Support, file = "Max Likelihood HGT Average.txt", sep = "     ")
write.table(Agreement.Support, file = "Combined HGT Average.txt", sep = "     ")
