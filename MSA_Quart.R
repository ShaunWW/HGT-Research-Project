
# Use this  section to load 'msa' package

# source("https://bioconductor.org/biocLite.R")
# biocLite("msa")

  # Examples
  # 
  # ## read sequences
  # filepath <- system.file("examples", "exampleAA.fasta", package="msa")
  # mySeqs <- readAAStringSet(filepath)
  # 
  # ## call unified interface msa() for default method (ClustalW) and
  # ## default parameters
  # msa(mySeqs)

library("msa", lib.loc="~/R/win-library/3.2")

#library("rphast", lib.loc="~/R/win-library/3.2")


Genes<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")
L <-length(Genes)



combo.vec<- matrix(data=0, nrow= 715, ncol = 4)
# ^ combo.vec hold the 13c4 combinations of geneID's and will be used later to read quartet files.
a<-1
b<-2
c<-3
d<-4

index<-1
# Generating Combo vec
while (a <= (L-3))
{
  while( b <= (L-2))
  {
    while(c <= (L-1))
    {
      
      while(d <= L)
      {
        combo.vec[index ,] <- c(Genes[a], Genes[b], Genes[c], Genes[d])
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

# WD with the quartets
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")


# Read Quarte table:
i<-1
while ( i<= length(combo.vec[,1]))
{
  Set<- combo.vec[i,]
  Quart.Table <- read.table(paste(c("quartGID", Set[1],"_",Set[2],"_", Set[3],"_",Set[4]), sep = "", collapse = ""), colClasses = "character")
  # ****^ Name may need to be changed, not verbose
  
  n.quart<- length(Quart.Table[[1]])
  
  j <- 1
  while( j <= n.quart ) # For each quart: get 4 seqs -> run MSA -> write MSA
  {
    seqs <- c( Quart.Table[[1]][[j]], Quart.Table[[2]][[j]], Quart.Table[[3]][[j]],Quart.Table[[4]][[j]])
    # ^ This should be the 4 sequences in the quartet
    
    # # test set
    # g1 <- "HWRTMMHPIAMRAAVPIPYSSAPIIAAMTMSRPVRRPPSVRSVTRSRRLFIASTWCASVRPISHGRPAYLMDVPGDAPVP"
    # g2 <- "MVVSQEPYTVRTALKKQRVLIFFPVWLVKSGVRWYVRTASQGANVLFSTFEGQTQERHAYGMCVYHIAHPGTETFLKT"
    # g3 <- "MRNIAKALTTNKPPEGARTLDPVLEGLRGFAARSPLRAQHSLSPPLPVAPLSLPSQPPHFPDASTRHTHHAEDEPSAGFG"
    # seqs<- c(g1, g2, g3)
      
    MSA <- msa::msa( seqs, method = "ClustalW", type = "protein")
    
    # FIXME
    # Need to write or save each MSA for jth quartet
    
    j <- j+1 #next quartet
  }
  
  i<- i+1 # next set
}


file:///C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits/orthologs_8_or_better/orthologs_1.txt