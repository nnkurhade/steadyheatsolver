library("plot3D")
library("plot3Drgl")
rm(list = ls())
phi <- read.table("F:/Studies/Mtech/Dept of Aerospace Engg/Sem2/ME6151/Project/readmesh/phi.txt", quote="\"", comment.char="")
x <- phi$V1
y <- phi$V2
z <- phi$V3
Z = matrix(0,length(x),length(y))
k<-1
for (i in 1:length(x)){
  for(j in 1:length(y)){
    if(i==j){
      Z[i,j]<-z[k]
      k<-k+1
    }
  }
}
scatter3D(x,y,z,type = "h")
plotrgl()