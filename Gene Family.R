# Gene Family Identifyer

library("micropan", lib.loc="~/R/win-library/3.2")
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")
Gene.Num<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")

###
# adjust hit.min for as desired 
hit.min <- 10

#for (a in 1:length(Gene.Num))
#{
  base.genome<- Gene.Num[2]
  gene.list <- c("1", "2", "3")
  index<-1
  # a<-a+1
  a<-1
  
  b<-2 ##
  for(b in 1:length(Gene.Num)) # This for loop generates table of entry Genomes's best matches in all genomes
    {   #For each comparison, read in genes from gene.num 'a'
    if(a != b)
      {
      base <- read.table(paste(c("Best Match_", base.genome, "&", Gene.Num[b],".txt"), sep="", collapse = ""), colClasses = "character")
      gene.list[index]<- base[1]
      index<- index+1 # should have N-1 lists of genes from 'a'
      }
    }
  
  
  # gene.list -> 12 lists of Best genes in genome 'a'
  
  #min.vec<-c("Sequence", "length.vec")
  max.vec<-c("Maximum", "Sequence")
  c <- 1
  while (c < length(Gene.Num))
    {
    #.vec[c]<- gene.list[[c]][1]
    max.vec[c]<- gene.list[[c]][length(gene.list[[c]])]
    
    c <- c+1
    }
  #alpha <- min(min.vec)
  beta  <- max(max.vec)
  

  seq <- 1
  Q<- "query"
  q.list <- c("store", "query")
  hit.list<- c("store", "hits")
  index<-1
  
  while (Q != beta) # while loop generates quereys spanning alpha and beta
    {
    Q<-  paste(c("GID", base.genome, "_seq", seq), sep="", collapse = "")
    hits<-0 # reset hits and index for Q
  
    #for each querey (Q) look for it in gene.list
    list <-1
    mark<-1

    while(list < length(Gene.Num)) # for each genome
      {
      # for ONE of twelve genomes
      # look for match
     
      j<- 1
     #if (seq < 10000)
        #{
        while( j < length(gene.list[[list]]))
          {#^ scrolls first 999 posibilities
          if (gene.list[[list]][j] == Q)
            {
            hits<- hits+1
            j <- length(gene.list[[list]]) # break while loop
            } 
          
          j<- j+1
          }
        #}
      
      # else if ( (seq > 999) && (seq < 10000))
      #   {
      #   for(j in 1000:9999)
      #     {#^ scrolls first 999 posibilities
      #       if (gene.list[[list]][j] == Q)
      #       {
      # 
      #         # length(gene.list[[list]])
      #         hits<- hits+1
      #         mark <- mark+1
      #       }
      # 
      #       j<- j+1
      #     }
      #   }
      # else
      #   {
      #   break
      #   }

      # else {j < 100000}
      
      list <- list+1 
      }
    
    # after each genome checked for Q
    if (hits >= hit.min) #threshold reaches, gene in family
      {
      hit.list[index]<- hits
      q.list[index]<- Q
      index<- index+1
      }
    
    seq <- seq+1
    }
  Gene.Fam <- data.frame(q.list, hit.list)
  
##___________________________________________________________________________________________-  
  # Script to prepare gene family for clust
    # need fasta format of genes in family
  
  library("seqinr", lib.loc="~/R/win-library/3.2")
  
  # Read in sequence file
  setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped")
  Chon<- "Chondrus_genbank_GID0001"
  genome<-read.fasta(file = Chon, seqtype = "AA", as.string = TRUE)
  
  # ## FIXME ***need to read data from file ***
  # setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")
  # family<- read.table()

  fam.seq<- c("Gene","family","sequences")
  
  #attr(genome[[1]], which = "name")
  #^ gets name of fasta sequence
  
  # Gene Family Generator    
  #   need to compare name to gene.fame[1]
  j<-1
  fam.index<- 1
  while ( j < length(genome) )
    {
    if(attr(genome[[j]], which = "name") == q.list[fam.index] )
      {
      fam.seq[fam.index]<- genome[[j]]
      fam.index<- fam.index+1
      }
    j<- j+1
    }
  Gene.Fam<-data.frame(q.list, fam.seq)
  
  for(i in 1:length(Gene.Fam[[1]]))
    {
    write.fasta(sequences = fam.seq[[i]], names = q.list[[i]], file.out = paste(c("GeneFam_", base.genome),sep="", collapse = ""),open= "a", as.string = FALSE )
    }
  
  
#_________________________________________________________________________________________________________________-  
  x<-c(strsplit(BLAST.vector[[1]], "_"))[[1]]
  substr(x[[1]],4, nchar(x[[1]]) )