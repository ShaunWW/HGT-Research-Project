Gene.Num<-c("0001", "0011", "0312", "0313", "1111", "1172", "1353", "1354", "2222", "3333", "4444", "5706", "06840")
L<- length(Gene.Num)

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