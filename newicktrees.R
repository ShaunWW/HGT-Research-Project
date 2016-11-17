# BACKGROUND: This program was used in a project where we were trying to determine 
# the evolutionary history of 14 different species of algae. This was done by 
# bootstrapping 100 "individuals" of each species from thier sequenced genome and 
# reconstructing the most likely tree.

# This program generates a matrix of 13 choose 4 combinations for every combination
# species, and reads in the coresponding bootstraped trees. It then figures out the
# most probably evolutionary history for the set of four, aka "the cherry". 



setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_trees")
# Generate matrix of possible 14c4 combos
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
#_________________________________________________________________________________________________________________________

# for each bootstrap file

# START P WHEREVER
p<-8
while ( P < 1970)
{
  # Read in bootstrapped trees
  tree.table <- read.table(file = paste(c("newicktrees_", P ,".txt"), sep ="", collapse = ""), stringsAsFactors = FALSE, sep = ";")
  trees<- tree.table[1]
  rm (tree.table)
  
  # Create and fill in tree data
  tree.data <- matrix( data = 0, nrow = 1001, ncol = 4) #colnames(c("ABCD", "AB_match", "AC_match", "AD_match"))
  
  # Check 1001 combinations for embedded quartetes
  for (i in 1:1001)
  {
    #read in quartet from combo.vec
    quart<-c(combo.vec[i,1],combo.vec[i,2], combo.vec[i,3], combo.vec[i,4])
    tree.data[i,1] <- paste(quart, sep = "", collapse = "_")
    check <- c(0,0)
    
    AB.match <- 0
    AC.match <- 0
    AD.match <- 0
    #read split for each tree (100)
    for (j in 1:100)
    {
      # Find cherry in Tree
      set.1 <- unlist(strsplit(trees[j,1], split = ""))
      #set.1 <- unlist(strsplit("((((1,3),(2,(7,5))),4),6);", split = ""))
      set.2 <- c("list", "delete", "( )")
      set.2.index <- 1
      t <- length(set.1)
      L <- 1
      while (L < t) # remove elements not part of quartet
      {
       if( set.1[L]== "(")
       { 
         # # If the last entry was NOT a "(", index
         # if ( set.2[set.2.index] != "("  && set.2.index != 1)
         # { set.2.index <- set.2.index +1}
         # # If the last entry WAS a "(",  do NOT index
        
         set.2[set.2.index]<- set.1[L]
         set.2.index <- set.2.index +1
       } else if (set.1[L] == ")")
        {
          # # If the last entry was NOT a ")", index
          # if ( set.2[set.2.index] != ")")
          # { set.2.index <- set.2.index +1}
          # 
          set.2[set.2.index]<- set.1[L] 
          set.2.index <- set.2.index +1
        } else if ( set.1[L]== ",")
        {
          set.2[set.2.index]<- set.1[L]
          set.2.index <- set.2.index +1
          
        } else # set.1[L] is a number
          {
          # check for two-digit number -> if next entry is a number
          if (set.1[L+1] != "(" && set.1[L+1] != ")" && set.1[L+1] != ",")
          { 
            # Executes if 2 digit number
            x <- as.numeric(paste(c(set.1[L], set.1[L+1]), sep ="", collapse = ""))
           
             L<- L+1 # Extra L index
          } else # One digit number 
            {
              x <- as.numeric(set.1[L])
            }
          # check if the number is in the given quartet
          if (is.element(x, quart))
          {
            set.2[set.2.index]<- x
            set.2.index<- set.2.index +1
          }
        }
        L <- L +1
      }
      set.2[set.2.index]<- ")"
      paste(set.2,sep = "", collapse = "")
      set <- strsplit(set.2, split = c(","))
      rm(set.1, set.2, set.2.index, L)
      
      # Find the parenthesis with exactly two numbers
      #   Those two numbers are the cherry
      g <-1
      start<-g
      end <- FALSE
      cherry <- c(0,0)
      while(end == FALSE && start < 1000)
      {
        start <- g
        num.count <- 0 # track the numbers in each parenthesis
        num.set <- c(0,0) # record these numbers
        #find open
        while( set[start] != "(")
        {
          start <- start +1
        }
        # find start, open
        open <-1 
        start <- start +1
        
        #find close
        while( open > 0)
        {
          if (set[start] == "(")
          {
            open <- open +1
          }
          
          if (is.element(set[start], quart))
          {
            num.count <- num.count +1  
            num.set[num.count] <- set[start]
          }
          
          if (set[start] == ")")
          {
            open <- open -1
          }
          
          start <- start +1
        }
        # When open == 0, the starting parenthesis is closed
        
        if ( num.count == 2) # the numbers contained are the cherry
        {
          end <- TRUE
          cherry <- num.set
        }
        
        g <- g+1
      }
     
      #check split and assign to correct bin
      check<- sort(as.numeric(cherry))
      
      if( check[1] == quart[1]  && check[2] == quart[2]) 
      { 
        AB.match = AB.match +1 
      } else if( check[1] == quart[3]  && check[2]== quart[4]) 
      { 
        AB.match = AB.match +1 
      } else if( check[1] == quart[1]  && check[2]== quart[3]) 
      { 
        AC.match = AC.match +1 
      } else if( check[1] == quart[2]  && check[2]==quart[4]) 
      {
        AC.match = AC.match +1 
      } else if( check[1] == quart[1]  && check[2]==quart[4]) 
      {
        AD.match = AD.match +1 
      } else if( check[1] == quart[2]  && check[2]==quart[3]) 
      {
        AD.match = AD.match +1
      } 
      
    }
    
    # write split data for last quart
    tree.data[i,2]<- AB.match
    tree.data[i,3]<- AC.match
    tree.data[i,4]<- AD.match
    
  }
  
  # Write tree.data to file
  write.table(tree.data, file = paste(c("TreeData_", P ,".txt"), sep ="", collapse = ""), sep = "   ")
  P <- P +1
}
