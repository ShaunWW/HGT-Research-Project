setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes")

# Name FASTA File to be read 
gald_sulph <- "Galdieria_sulphuraria_GID111.fsa"
font <-"Fonticula_alba_GID123.fsa.fsa"
tryp <-"Trypanosoma_brucei_GID789.fsa"
para <-"Cryptomonas_paramecium_GID456.fsa.fsa"

# Load seqinr
library("seqinr", lib.loc="~/R/win-library/3.2")


# Read in FASTA files
fasta.res.GALD <- read.alignment(file = gald_sulph,format = "fasta")
fasta.res.FONT <- read.alignment(file = font,format = "fasta")
fasta.res.TRYP <- read.alignment(file = tryp,format = "fasta")
fasta.res.PARA <- read.alignment(file = para,format = "fasta")

# Translate FASTA files REQUIRE SEQUIN TO BE DISBALED to translate
detach("package:seqinr", unload=TRUE) # disable seqinr

para_PROT <- vector(mode="character", length = length(fasta.res.PARA[[2]]))
for(i in 1:length(fasta.res.PARA[[2]])){
  para_PROT[[i]]<-translate(fasta.res.PARA[[3]][[i]])
}

# Turn on seqinr
library("seqinr", lib.loc="~/R/win-library/3.2")


# Write Protien to File
# first create file to append -> "Gald_Sulph_PROT_GID111 -> MUST have GID names
for(i in 1:length(fasta.res.PARA[[2]])){
  write.fasta(sequences = para_PROT[[i]], names = fasta.res.PARA[[2]][[i]], file.out = "para_PROT_GID456",open= "a", as.string = FALSE )
}

# Define protien file variable
FONT <-"font_PROT_GID123"
GALD <-"Gald_Sulph_PROT_GID111"
TRYP <-"tryp_PROT_GID789"
PARA <-"para_PROT_GID456"


# BLAST Protien Files
blastAllAll(in.files = c(PARA,FONT,GALD,TRYP), out.folder = "BLAST", e.value=0.0001, job = 1, verbose = TRUE)

setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/BLAST")


filenames <- c("GID789_vs_GID789.txt", "GID789_vs_GID456.txt","GID789_vs_GID123.txt","GID789_vs_GID111.txt")

blast.distance <-bDist(filenames, e.value = 0.0001)

x<-readBlastTable("GID111_vs_GID111.txt")



# Prepping Protien Genome Files:
setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped")
library("micropan", lib.loc="~/R/win-library/3.2")

# Define file varibales
GaldSul<-"Protein .fasta"

# Prep files for BLAST
panPrep(in.file = GaldSul, GID.tag = "GID1994",out.file = "Gald_Sulph", protein = TRUE)




  #Naming File variables
Timsp <-"Timspurckia_oligo_GID1172"
Rhodo <-"Rhodosorus_marinus_GID0011"
Porphy<-"Porphyridium_aerugineum_GID0313"
Porphy2<-"Porphyra_umbilicalis_GID06840"
Erymad<-"Erythrolobus_madagascarensis_GID1354"
Eryaus<-"Erythrolobus_australicus_GID1353"
Compso<-"Compsopogon_coeruleus_GID0312"
Chon<- "Chondrus_genbank_GID0001"
Calli<- "Calliarthron_tuberculosum_GID4444"
Cyanid<- "cyanidioschyzon_GID3333.organelle"
Gald_Phl<-"Galdieria_phlegrea_GID2222"
Pryo<-"Pyropia_haitanensis_GID5706"
Porphy3<-"Porphyridium_purpureum_GID1111"
GaldSul<-"Gald_Sulph_GID1994"

BLAST_Files <-c(GaldSul,Timsp,Rhodo,Porphy,Porphy2,Eryaus,Erymad,Compso,Chon, Calli, Cyanid,Gald_Phl,Pryo,Porphy3)

blastAllAll(in.files = BLAST_Files,out.folder = "red_algalBLAST",e.value =0.0001,job=1,verbose=TRUE)

setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/Prepped/red_algalBLAST")

bDist_File<-c("GID1354_vs_GID1354.txt", "GID1354_vs_GID1353.txt","GID1354_vs_GID1172.txt","GID1354_vs_GID06840.txt","GID1353_vs_GID1353.txt","GID1353_vs_GID1172.txt","GID1353_vs_GID1354.txt","GID1353_vs_GID06840.txt","GID06840_vs_GID1172.txt","GID06840_vs_GID06840.txt","GID06840_vs_GID1353.txt","GID06840_vs_GID1354.txt","GID1172_vs_GID1172.txt","GID1172_vs_GID1353.txt","GID1172_vs_GID1354.txt","GID1172_vs_GID06840.txt")

quart1234<-bDist(blast.files = bDist_File, e.value = 0.0001, verbose = TRUE)
