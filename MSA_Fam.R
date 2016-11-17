
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


#WD needs to have ppanpreped sequence files in gene families.

library("msa", lib.loc="~/R/win-library/3.2")
library("micropan", lib.loc="~/R/win-library/3.2")
library("seqinr", lib.loc="~/R/win-library/3.2")
library("Rphylip", lib.loc="~/R/win-library/3.2")
library("phangorn", lib.loc="~/R/win-library/3.2")
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits/formatted_orthologs")
Phylip.Path <- file.path("C:", "Users", "Shaun", "Desktop", "HGT Research Project", "phylip-3.695", "exe" )
#setwd(Phylip.Path)
#Genes<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")
#L <-length(Genes)

#1970 families of +8
I<-981
while ( I <= 1970 )
{
  
  file<- paste(c("orthologs_", I, ".txt"), sep="", collapse = "")
  fam.table<- read.table(file = file, header = FALSE , col.names = c("ID", "Sequence"), colClasses = "character")
  
  seqs<- fam.table[,2]
  names<- fam.table[,1]
  
  # Create empty matrix
  leng<- c(0,0,0)
  for (i in 1:length(names))
  {
    leng[i]<- nchar(seqs[i])
  }
  dim <- max(leng)
  mat.table <- matrix( data = NA, ncol = dim, nrow = length(seqs) )
  
  # Fill matrix with chars
  for (i in 1:length(seqs))
  {
    chars <- strsplit( seqs[i], fixed = TRUE, split = "")
    for( j in 1:length(chars[[1]]))
    {
      mat.table[i,j]<- chars[[1]][j] 
    }
  }
  dimnames(mat.table) <- list(names)
  # Matrix Done
  rm(i,j,names, seqs, dim, leng, chars)
  
  # Write phyDay file for ortholog family
  phyloData<-phyDat(data = mat.table, type = "AA")
  #write.phyDat(phyloData, file = paste(c("phyDat_", I, ".txt"), sep="", collapse = ""), format = "phylip")
  
  # Convert phyDat into proSeq
  proteinInput <- as.proseq(phyloData)
  
  # Run seqboot
  boot.data<-Rseqboot( proteinInput, method = "bootstrap")
  file.rename("outfile", paste(c("BootData_", I, ".txt"), sep="", collapse = ""))

  I <- I+1
}



{
  # Need to write in Formated fil
  
  # #Shorten Names, Get longest
  # file.create(paste(c("PhydatForm_", I, ".txt"), sep="", collapse = ""))
  # for (i in 1:length(seqs))
  # {
  #   x <- strsplit(names[i], "_")[[1]]  
  #   names[i]<- x[1]
  #   
  #   leng[i]<- nchar(seqs[i])
  #   format(names[i], width = 10)
  #   
  #   write(c(format(names[i], width = 10), seqs[i]), file = paste(c("PhydatForm_", I, ".txt"), sep="", collapse = ""), append = TRUE , sep = "")
  #   
  # }
  # max(leng)
  # file<-paste(c("PhydatForm_", I, ".txt"), sep="", collapse = "")
  # #dat<- read.phyDat(file, type = "AA")
  # #attr(MSA, which = "width")
  # #MSA@width
  # #msa::.__C__MsaMetaData
  # 
  #write.f
}

 g1 <- "HWRTMMHPIAMRAAVPIPYSSAPIIAAMTMSRPVRRPPSVRSVTRSRRLFIASTWCASVRPISHGRPAYLMDVPGDAPVP"
 g2 <- "MVVSQEPYTVRTALKKQRVLIFFPVWLVKSGVRWYVRTASQGANVLFSTFEGQTQERHAYGMCVYHIAHPGTETFLKT"
 g3 <- "MRNIAKALTTNKPPEGARTLDPVLEGLRGFAARSPLRAQHSLSPPLPVAPLSLPSQPPHFPDASTRHTHHAEDEPSAGFG"
 seqs<- c(g1, g2, g3)







combo.vec<- matrix(data=0, nrow= 715, ncol = 4)
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









# # Read Quarte table:
# i<-1
# while ( i<= length(combo.vec[,1]))
# {
#   
#   
#   Set<- combo.vec[i,]
#   Quart.Table <- read.table(paste(c("quartGID", Set[1],"_",Set[2],"_", Set[3],"_",Set[4]), sep = "", collapse = ""), colClasses = "character")
#   # ****^ Name may need to be changed
#   
#   #Quart/OP <- # i'th line of table
#     
#   i<- i+1
# }
# 
# 
