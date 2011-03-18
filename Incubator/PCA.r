mydata=read.table("ped.txt", header=T)
UIdx=which(mydata[,3]==0)
RIdx=which(mydata[,3]!=0)
d=mydata[,-c(1:6)]
D=scale(d, T, T)

UD=D[UIdx,]
pca=prcomp(UD)
plot(pca$x[,1], pca$x[,2])

E=eigen(cor(UD))
pc1=D%*%E$vector[,1]
pc2=D%*%E$vector[,2]

plot(pc1,pc2)

myphe=read.table("phe.txt", header=T)
model=glm(family=binomial, myphe[,3]~pc1)
