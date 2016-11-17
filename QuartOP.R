# This program takes the list of four geneomes and finds QuartOP's


library("micropan", lib.loc="~/R/win-library/3.2")

#   *** WORKING DIRECTORY NEEDS TO CONTAIN "Best Match_" FOR EACCH PAIR OF GENOMES ***
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")


# GENERATE LISTS
# Pick four Genomes for the Quartet


# Genes corresponds to the number associated with each geneome, they are indexed 1-13, and are referenced below
Genes<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")
L <-length(Genes)


# These four genes will be the first document made by the 1st iteration. 
# 
# #* In this case it is set to "0001", 1354, "4444", and "5706"
# gene1<-1
# gene2<-8
# gene3<-11
# gene4<-12

# I have set this one where my documents leave off. 
gene1<-2
gene2<-3
gene3<-7
gene4<-9


# The outer section is the algorithm I wrote to choose 4 from 13
while(gene1 <= (L-3))
{
  while(gene2 <= (L-2))
  {
    while(gene3 <= (L-1))
    {
      while(gene4 <= length(Genes))
      {#   *** This Section is the Program Core ***
        
        Gene.Num<-c(Genes[gene1], Genes[gene2], Genes[gene3], Genes[gene4])
        base.genome <- Gene.Num[1] # set base genome
        gene.list <- c("1", "2", "3")
        index<-1
        
        for(b in 2:length(Gene.Num)) # This for loop generates table of entry Genomes's best matches in all genomes
        {   
          #For each comparison, read in genes from gene.num 'a'
          
          # $$$ - base is read from a table by the name generator below, these are the "Best Match_" files I sent
          base <- read.table(paste(c("Best Match_", base.genome, "&", Gene.Num[b],".txt"), sep="", collapse = ""), colClasses = "character")
          
          
          
          gene.list[index]<- base[1]
          index<- index+1 # should have N-1 lists of genes from 'a'
        }
        # clean up useless vars
        rm(index,b)
        
        #_______________________________________________________________________________________________#
        # STEP #1: generate list of all genes in base genome that appear in AB,AC,AD
        
        #find beta - maximal gene in A
        length.vec<- c( gene.list[[1]][length(gene.list[[1]])],gene.list[[2]][length(gene.list[[2]])],gene.list[[3]][length(gene.list[[3]])] )
        beta<-max(length.vec)
        rm(length.vec)
        
        seq <- 1
        Q<- "query"
        q.list <- c("store", "query")
        index<-1
        
        while (Q != beta) # while loop generates quereys spanning alpha and beta
        {
          Q<-  paste(c("GID", base.genome, "_seq", seq), sep="", collapse = "")
          hits<-0 # reset hits and index for Q
          
          #for each querey (Q) look for it in gene.list
          list <-1
          
          while(list < length(Gene.Num)) # for each genome
          {
            # for ONE of 3 genomes
            # look for match
            j<- 1
            while( (gene.list[[list]][j] != Q) && (j < length(gene.list[[list]])) )
            {
              j<- j+1
              #After this loop, EITHER j=length(list) or Q is a hit
            }
            
            if(gene.list[[list]][j] == Q)
            {
              hits<- hits+1    
            }
            list <- list+1 
          }
          
          # after each genome checked for Q
          if (hits == 3) # gene found in all 3 lists
          {
            q.list[index]<- Q
            index<- index+1
          }
          
          seq <- seq+1
        }
        
        #q.list Hold list from step 1
        rm(seq, Q, list, j, index, hits)
        #________________________________________________________________________________________________________#
        #STEP #2: Find B_x, C_y, and D_z for each gene in A
        
        mm <- matrix(0, length(q.list), 4)
        mm[,1] <- q.list
        b<-2
        for(b in 2:length(Gene.Num)) # look for all matches in one list at a time
        {
          #For each comparison, read in genes from gene.num 'a'
          
          # $$$ This is also reading in a "Best_match folder"
          base <- read.table(paste(c("Best Match_", base.genome, "&", Gene.Num[b],".txt"), sep="", collapse = ""), colClasses = "character")
          
          
          
          match.list<- c("matches"," For A")
          j<-1
          while (j <= length(q.list))
          {
            A_n<- q.list[j]
            index<-1  
            
            # Look for A_n in base
            f<-1
            while(base[[1]][f] != A_n)
            {
              f<-f+1
            }
            
            match.list[[j]]<- base[[2]][f]
            j<- j+1
          }
          mm[,b]<-match.list
        }
        
        #a.matches contains the matches in A for genes that matched in B,C,&D
        a.matches<- data.frame(mm)
        
        rm(base, f, j, index, b, A_n)
        
        #________________________________________________________________________________________________________#
        #STEP #3: Check A's matches against each other
        # Define mm[,1]= A, mm[,2]= B, mm[,3]= C, mm[,4]= D
        
        # Need BC,BD,CD tables
        table.BC <- read.table(paste(c("Best Match_", Gene.Num[2], "&", Gene.Num[3],".txt"), sep="", collapse = ""), colClasses = "character")
        table.BD <- read.table(paste(c("Best Match_", Gene.Num[2], "&", Gene.Num[4],".txt"), sep="", collapse = ""), colClasses = "character")
        table.CD <- read.table(paste(c("Best Match_", Gene.Num[3], "&", Gene.Num[4],".txt"), sep="", collapse = ""), colClasses = "character")
        
        
        Quart.list<- matrix(0, length(q.list), 4)
        index <-1
        j<-1
        while(j< length(mm[,1]))
        {
          Quart.OP<-FALSE
          # matches to A.n, check against each other
          B.x<-mm[j,2]
          C.y<-mm[j,3]
          D.z<-mm[j,4]
          
          # find B.x in BC 
          bc<-1
          while( (table.BC[[1]][bc] != B.x) && (bc < length(table.BC[[1]])) )
          { bc <- bc+1 }
          
          # find B.x in BD
          bd<-1
          while( (table.BD[[1]][bd] != B.x) && (bd < length(table.BD[[1]])) )
          { bd <- bd+1 }
          
          
          if( (table.BC[[2]][bc]==C.y)  && (table.BD[[2]][bd] == D.z) )
          {
            # If this is true, C.y and D.z are Bx's best match, check if C.y matches with D.z
            cd<-1
            while( (table.CD[[1]][cd] != C.y) && (cd < length(table.CD[[1]])) )
            { cd <- cd+1 }
            
            if(table.CD[[2]][cd]==D.z)
            {Quart.OP <- TRUE}
          }
          
          if (Quart.OP == TRUE)
          {
            # If this is true, A.n B.x C.y D.z are QuartOP!
            Quart.list[index,]<- mm[j,]
            index<- index+1
          }
          
          
          
          j<- j+1
        }
        
        write.table(Quart.list, file = paste(c("quartGID",Gene.Num[1],"_",Gene.Num[2],"_",Gene.Num[3],"_",Gene.Num[4]), sep = "",collapse = ""), sep = "     ", append = FALSE)
        
        
          # *** This is the End of the Core ***
        
        gene4<- gene4+1
      }
    
      gene3<- gene3+1
      gene4<- gene3+1
    }
    
    gene3<- gene2+1
    gene2<- gene2+1
  }
  
  gene2<- gene1+1
  gene1<- gene1+1
}

