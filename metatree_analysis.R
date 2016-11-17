# This program combines and read the 5 different meta trees, coresponign to the 
# level of support (50%, 60%...90%). 

setwd("C:/Users/Shaun/Desktop/HGT Research Project/Prepped Geneomes/red_algal_proteins.est-derived (1)/red_algal_proteins.est-derived/bootstrapped_splits")

best_split.mat <- matrix(nrow = 1001, ncol= 7)

for (i in 5:9)
{
  meta =       paste("Metafile_",i,"0.txt", sep = "", collapse = "")
  best_split = paste("Quartet Best Split_",i,"0.txt", sep = "", collapse = "")
  
  meta.table =       read.table( file = meta, stringsAsFactors = FALSE)
  best_split.table = read.table( file = best_split, stringsAsFactors = FALSE)
  
  for ( j in 1:1001)
  {
    best_split.mat[j,i-3] = best_split.table[j,2]
    best_split.mat[j,1] =   best_split.table[j,1]
  }
}

for ( j in 1:1001)
{
  # check if all entries are the same
  
}

write.table(best_split.mat , file = "Best Split Matrix.txt")
