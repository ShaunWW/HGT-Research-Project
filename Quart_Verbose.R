# THis Program reads in a Quart list and writes a verbose file.


# Prepping Protien Genome Files:
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped")
library("seqinr", lib.loc="~/R/win-library/3.2")


#Naming File variables
T.0001<- "Chondrus_genbank_GID0001"
Rhodo <-"Rhodosorus_marinus_GID0011"
Compso<-"Compsopogon_coeruleus_GID0312"
Porphy<-"Porphyridium_aerugineum_GID0313"
Porphy3<-"Porphyridium_purpureum_GID1111"
Timsp <-"Timspurckia_oligo_GID1172"
Eryaus<-"Erythrolobus_australicus_GID1353"
Erymad<-"Erythrolobus_madagascarensis_GID1354"
Gald_Phl<-"Galdieria_phlegrea_GID2222"
Cyanid<- "cyanidioschyzon_GID3333.organelle"
Calli<- "Calliarthron_tuberculosum_GID4444"
Pryo<-"Pyropia_haitanensis_GID5706"
Porphy2<-"Porphyra_umbilicalis_GID06840"

T.0001<- read.fasta(file = Chon, seqtype = "AA", as.string = FALSE, seqonly = FALSE)
T.0011<- read.fasta(file = Rhodo, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.0312<- read.fasta(file = Compso, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.0313<- read.fasta(file = Porphy, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.1111<- read.fasta(file = Porphy3, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.1172<- read.fasta(file = Timsp, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.1353<- read.fasta(file = Eryaus, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.1354<-read.fasta(file = Erymad, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.2222<-read.fasta(file = Gald_Phl, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.3333<-read.fasta(file = Cyanid, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.4444<-read.fasta(file = Calli, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.5706<-read.fasta(file = Pryo, seqtype = "AA", as.string = TRUE, seqonly = FALSE)
T.06840<- read.fasta(file = Porphy2, seqtype = "AA", as.string = TRUE, seqonly = FALSE)

rm(Chon, Rhodo, Compso, Porphy, Porphy3, Timsp, Eryaus, Erymad, Gald_Phl, Cyanid, Calli, Pryo, Porphy2)

Genes<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")
L <-length(Genes)

# setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Unprepped")
# i<-1
# while(i <= L)
# {
#   tab<- read.alignment( file= paste(c("GID",Gene.Num[i],".txt"), sep = "",collapse = ""), format = "fasta") 
# }

# Need a martix of lists
#Protien.file.ID<- matrix(data = NA, nrow = 19167, ncol = 13)


setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/BLAST Top Hits")
gene1<-1
gene2<-2
gene3<-3
gene4<-4

# Outer Loop -> 13 choose 4
while(gene1 <= (L-3))
{
  while(gene2 <= (L-2))
  {
    while(gene3 <= (L-1))
    {
      while(gene4 <= length(Genes))
      {
        # MAIN LOOP BEGINS ________________________________________________________________________________________
        
        # Generate necessary files
        Gene.Num<-c(Genes[gene1], Genes[gene2], Genes[gene3], Genes[gene4])
        quart.file<- read.table(file = paste(c("quartGID",Gene.Num[1],"_",Gene.Num[2],"_",Gene.Num[3],"_",Gene.Num[4]), sep = "",collapse = ""), as.is = TRUE)
        
        # While loop to find number of quartes, make matrix
        j<-1
        while (quart.file[[1]][j] != 0)
        {
          j<- j+1  
        }
        num.quarts<- j-1
        
        Verb.Quart<- matrix( data = "DNA_Sequence", nrow = num.quarts, ncol = 4)
        
        #paste(c("quartGID",Gene.Num[1],"_",Gene.Num[2],"_",Gene.Num[3],"_",Gene.Num[4]), sep = "",collapse = "")
        
        i<- 1
        library("seqinr", lib.loc="~/R/win-library/3.2")
        while(i <= length(Gene.Num))
        { # This while loop runs for each of four genomes
        
          Table<- read.fasta(paste(c("T_",Gene.Num[i]), sep = "",collapse = "") , seqtype = "AA", as.string = TRUE, seqonly=FALSE)
          j<-1
          while (j <= num.quarts)
          { 
            #This while loop runs through all entries of each genome list
            Q<- quart.file[[i]][j]
            x<-1
            while(attr(Table[[x]], which = "name") != Q )
            {
              x<-x+1
              # When this loop ends, x stores index in table with sequence
            }
            # need to translate and write sequence
            
            detach("package:seqinr", unload=TRUE) # disable seqinr
            library("micropan", lib.loc="~/R/win-library/3.2")
            
            Verb.Quart[j,i]<-translate( seq=Table[[x]][1] , sens = "R")
            
            j<- j+1
          }
          
          i<-i+1
        }
        
        # After this loop, we have our sequences, need to translate
        
        
        gene4<- gene4+1
        # MAIN LOOP ENDS_________________________________________________________________________________________
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

rm(gene1, gene2, gene3, gene4)
