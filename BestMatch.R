# Quartet: Eukaryotic Algae 
#   HGT know to be prevelant
#   Studying the RATE of HGT

# given genomes 1,2,3,4 containing genes A,B,C,D, if
# 
# the best match for 1A in genome 2 is 2B  ((GOT IT))
# the best match for 1A in genome 3 is 3C
# the best match for 1A in genome 4 is 4D
# the best match for 2B in genome 3 is 3C
# the best match for 2B in genome 4 is 4D
# the best match for 3C in genome 4 is 4D
# 
# and all 6 reverse statements hold as well, then we call (1A,2B,3C,4D) a QuartOp.


# Six files for maximal gene matches for all 6 genome pairs

# quartet with high bootstrap support conflicts with consense -> MIGHT BE HGT
# SCRIPT TO FIND Mutual Maximum Matches (MMM) for two genomes

library("micropan", lib.loc="~/R/win-library/3.2")
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/red_algalBLAST")

# List of BLAST files to find best match
BLAST.vector <- c("GID1172_vs_GID1172.txt","GID0011_vs_GID0011.txt","GID0313_vs_GID.txt","GID06840_vs_GID1994.txt","GID1353_vs_GID1994.txt", "GID1354_vs_GID1994.txt", "GID0312_vs_GID1994.txt", "GID0001_vs_GID1994.txt", "GID4444_vs_GID1994.txt", "GID3333_vs_GID1994.txt", "GID5706_vs_GID1994.txt", "GID2222_vs_GID1994.txt", "GID1111_vs_GID1994.txt", "GID####_vs_GID1994.txt")

for (r in 1:length(BLAST.vector))
{
  AB<-readBlastTable(BLAST.vector[[r]])
  
  x<-c(strsplit(BLAST.vector[[r]], "_"))[[1]]
  genome1 <- x[[1]]
  genome2 <- x[[3]]

  align <- 0 # setting initial alignment value
  seq<-1 # start at the first sequence
  best_match_row <-0
  
  # define data table vectors
  Query<- c("Sequence","query", "yolo")
  Best.Match <- c("Best", "Match","y0")
  Percent.match<- c(0, 100, 69)
  List_entry<-1
  i = 1;
  
  while (i <= length(AB[[1]]))
    {
    # Code to find the maximal alignment with seqA_j
    if(AB[[i,1]]== paste(c(genome1,"_seq",seq), sep="", collapse = ""))
    {
      #^ If this is true, AB[[i,3]] is the alignment of seqA_j with SeqB_x
      
      # Need to test this agreement with other matches of seqA_j
      if (AB[[i,3]] > align)
      {
        align <- AB[[i,3]]
        best_match_row <- i # stores the row of the best reading
      }
      i = i+1
    }
    # seq_j readings done,
    else if (best_match_row > 0)  
    {    
      # Save Quesry, Best hit, and %ID
      Query[[List_entry]]<- AB[[best_match_row,1]]
      Best.Match[[List_entry]]<- AB[[best_match_row,2]]
      Percent.match[[List_entry]] <- align
        
    
      List_entry <- List_entry+1     
      align <- 0
      best_match_row<-0
    }
    else 
    {
      # If this is true, No queries found for seq_j
      seq<- seq+1
    }
  }
  Query[[List_entry]]<- AB[[best_match_row,1]]
  Best.Match[[List_entry]]<- AB[[best_match_row,2]]
  Percent.match[[List_entry]] <- align
  
  
  MMM_AB <- data.frame(Query,Best.Match,Percent.match)
  
  write.table(MMM_AB, file = paste(c("MMM_",genome1,"_vs_", genome2), sep="", collapse = "") , quote = FALSE,  sep = "       ", col.names = TRUE)

}
#---------------------------------------------------------------------------------------

# Compare the two tables to find MMM's between the two
#   Longer genome goes first (AB)

library("micropan", lib.loc="~/R/win-library/3.2")
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")
Gene.Num<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")

for (a in 7:length(Gene.Num))
  {
  base.genome<- Gene.Num[a]
  a<-a+1
  
  for(b in 1:length(Gene.Num))
    {
    if(base.genome != Gene.Num[b])
      {
      MMM_AB <- read.table(paste(c("MMM_GID", base.genome, "_vs_GID", Gene.Num[b],".txt"), sep="", collapse = ""), colClasses = "character")
      MMM_BA <- read.table(paste(c("MMM_GID", Gene.Num[b], "_vs_GID", base.genome,".txt"), sep="", collapse = ""), colClasses = "character")
  
      seq <- 1
      j<- 1
      list.index <- 1
      MMM.list.AB<- c("hi", "dude")
      MMM.list.BA<- c("hi", "dude")
      
      
      while (j <=length(MMM_AB[[1]]))
      {
        if(MMM_AB[[j,1]]== paste(c("GID", base.genome, "_seq", seq), sep="", collapse = ""))
        {
          j.MMM<- FALSE
          # ^ If this is true, look for j's match in BA, check if its match is seq_j
          seq.j <- MMM_AB[[j,1]] # store value of seq_j's name
        
          t<-1
          while (t < length(MMM_BA[[1]]))
          {
            if(MMM_BA[[t,2]] == seq.j) # look for j's best match in BA
            {
              j.MMM<- TRUE # seq_j is MMM
              break
            }
            t<- t+1
          }
      
          
          if (j.MMM == TRUE) # j is MMM
          {
            MMM.list.AB[[list.index]]<- MMM_AB[[j,1]]
            MMM.list.BA[[list.index]]<- MMM_AB[[j,2]]
            
            list.index <- list.index+1
          }
          
          j<- j+1;
        }  
        else { seq <- seq+1} # no seq.j found, increment
      }
        
      MMM.list <- data.frame(MMM.list.AB, MMM.list.BA)
        
      write.table(MMM.list, file = paste(c("Best Match_", base.genome, "&", Gene.Num[b],".txt"), sep="", collapse = ""), quote = FALSE,  sep = "       ", col.names = TRUE)
      }
      
    b<- b+1
    }
  }

#________________________________________________________________________________

  