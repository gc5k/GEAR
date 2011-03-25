mydata=read.table("ped.txt", header=T)
UIdx=which(mydata[,3]==0)
RIdx=which(mydata[,3]!=0)
d=mydata[,-c(1:6)]
D=scale(d, T, T)

UD=D[UIdx,]
pca=prcomp(UD)
plot(pca$x[,1], pca$x[,2])

kc1=D%*%pca$rotation[,1]
kc2=D%*%pca$rotation[,2]

plot(kc1, kc2)

myphe=read.table("phe.txt", header=T)
model=glm(family=binomial, myphe[,3]~kc1)
plot(myphe[,3]-model$fitted.values)
