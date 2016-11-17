
# current directory should contain folders "bootstrapped_trees" (input), 
# and "bootstrapped_splits" (output, empty) 


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

## The following is a function. It doesn't run by itself, but I can call it, and it'll then run, and
## return a value.

## This whole script sort of reads backwards, because I have to define a function before I can use it.


## record and return all numbers in current tree
getleaflist <- function(){
leaflist=c();
i=1;
curmax=0;
while (i<length(temp)){
## if current spot is not part of a number, skip
if ((temp[i] == ")") | (temp[i] == "(") | (temp[i] == ",")){
i=i+1;
}
curstring=c();

## when we find a number, convert it to integer form, and store it
while ((temp[i] != ")") & (temp[i] != "(") & (temp[i] != ",")){
curstring=paste(curstring,temp[i],sep="");
i=i+1;
}
myx=strtoi(curstring);
leaflist=c(leaflist,myx);
}
return (leaflist);
}


## Record all the quartets we can build by pairing two leaves in current sack 
## (that is, two leaves contained within a particular pair of parentheses) and 
## two leaves not in sack. The valid split is always the one in which the two quartets
## in the sack form a cherry. This method, when applied globally, will record some 
### quartets more than once, but that's okay. 

processsack <- function(sack){
a=length(sack);
for (i in 1:(a-1)){
for (j in (i+1):a){
leftpair=c(sack[i],sack[j]);
for (k in 1:(nleaf-1)){
if (sum(1*(sack==leaflist[k]))==0){
for (ell in (k+1):nleaf){
if (sum(1*(sack==leaflist[ell]))==0){
rightpair=c(leaflist[k],leaflist[ell]);
fourset=c(rightpair,leftpair);
if (sum(1*(leftpair==min(fourset)))>0){
pairone=leftpair;
pairtwo=rightpair;
}
else {
pairone=rightpair;
pairtwo=leftpair;
}
if (max(pairone)<min(pairtwo)){
cursplit=1;
}
else if (max(pairone)<max(pairtwo)){
cursplit=2;
}
else {cursplit=3;}
inds=sort(c(leftpair,rightpair));
B[inds[1],inds[2],inds[3],inds[4],cursplit] <<- 1;
}}}}}}
}

## The function "explore" is the main function.  It examines the whole area between
## two parentheses, records all the leaves inside, and forms all quartets formed from
## two leaves inside the parentheses and two outside. 
## If it finds a pair of parentheses inside the pair it's considering, it calls itself 
## on that interioor pair. 
## NOTE: When a function calls itself, that's known as "recursion." Recursion is a really powerful
## frequently quite elegant technique.

explore <- function(myinput){

## The assignment operator <<- means, assign this to a GLOBAL variable. (By default, all variables
## in a function are local, which means they're totally distinct from variables run outside the
## function, even if the two variables have the same name.
## Global variables are considered terrible programming style in all languages. They can
## indeed get you into lots of trouble, but as long as one is careful, I consider that they're okay.
## I try not to use them, but do so when it's really convenient.
curspot <<- myinput+1;

## Gather up all the leaves between this pair of parentheses, and store them in a sack.
cursack=c();
sackcount=0;

## Keep going till we reach the ) matching the ( which we entered at.
while (temp[curspot] != ")"){
if (temp[curspot]=="("){
## If we encounter a pair of parentheses, get all the leaves between them, and store in sack
if (sackcount==0){
leftsack=explore(curspot);
if (length(leftsack)>1){
## Record all quartets formed from two elements in left-hand sack, and two outside it
processsack(leftsack);
}
sackcount=sackcount+1;}
else {
## If sackcount=1, that means we've already explored the left-hand branch. So store whatever's
## in the right-hand branch in a new sack.
rightsack=explore(curspot);
if (length(rightsack)>1){
processsack(rightsack);
}
}

}
## we emerge from the case where we encountered an interior pair of parentheses

## if we hit a comma, just move along
else if (temp[curspot]==","){
curspot <<- curspot+1;
}
else{
## If we don't hit , ( or ), we must be at a number. Store it.
curstring=temp[curspot];
curspot <<- curspot+1;
## Take care of two-digit numbers
while ((temp[curspot] != ")") & (temp[curspot] != "(") & (temp[curspot] != ",")){
curstring=paste(curstring,temp[curspot],sep="");
curspot <<- curspot+1;
}

## Convert number of integer
myx=strtoi(curstring);

## Store number in sack. In this case there's no need to compute quartets based on the sack,
## since it only contains one element.
if (sackcount==0){
leftsack=myx;
sackcount=sackcount+1;}
else {
rightsack=myx;
}
}

}

## Combine both sacks, advance past the closing ), return sack, done.
sack=c(leftsack,rightsack);
curspot <<- curspot+1;

return(sack);
}



### THE ACTUAL SCRIPT STARTS HERE - EVERYTHING ABOVE IS JUST FUNCTIONS, WHICH WILL
##  BE CALLED BELOW.
P<-1
while ( P < 1971)
{
  tree.data <- matrix( data = 0, nrow = 1001, ncol = 4) #colnames(c("ABCD", "AB_match", "AC_match", "AD_match"))

  tree.table <- read.table(file = paste(c("bootstrapped_trees/newicktrees_", P ,".txt"), sep ="", collapse = ""), stringsAsFactors = FALSE, sep = ";")
  trees<- tree.table[1]
  rm (tree.table)
  temp <<- unlist(strsplit(trees[1,1], split = ""));
  ## record all leaves in first tree
  leaflist=getleaflist();
  nleaf <<- length(leaflist);
  ## We'll record our splits for the current file in A. 
  A <<- array(0,dim=c(14,14,14,14,3));

    for (j in 1:100){
	quartcount <<- 0;
    temp <<- unlist(strsplit(trees[j,1], split = ""));
	curspot=1;
	quartcount=0;
	# We'll use the array B to record the splits in current tree
	B <<- array(0,dim=c(14,14,14,14,3));
	# The following line does all the work - it explores current tree and records the quartets
	explore(curspot);
	# Add current splits to split total in A
	A <<- A+B;



}
  for (i in 1:1001){
    quart<-c(combo.vec[i,1],combo.vec[i,2], combo.vec[i,3], combo.vec[i,4]);
    tree.data[i,1] <- paste(quart, sep = "", collapse = "_")
	
  for (t in 2:4){
  # for current quartet, record the number of times each split occurred
	tree.data[i,t] = A[combo.vec[i,1],combo.vec[i,2], combo.vec[i,3], combo.vec[i,4],t-1];
  }
  }
  # write to file
    write.table(tree.data, file = paste(c("bootstrapped_splits/TreeData_", P ,".txt"), sep ="", collapse = ""), sep = "   ")
P <- P+1;
} 


