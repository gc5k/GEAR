mydata=read.table("F_ped.txt", header=T)
UIdx=which(mydata[,3]==0)
RIdx=which(mydata[,3]!=0)
d=mydata[,-c(1:6)]
D=scale(d, T, T)

UD=D[UIdx,]
pca=prcomp(UD)
plot(pca$x[,1], pca$x[,2])

kcc=read.table("kcc_ped.txt", header=T)
kc=kcc[,-c(1:6)]
kc=scale(kc, T, T)
kc1=kc%*%pca$rotation[,1]
kc2=kc%*%pca$rotation[,2]

plot(kc1, kc2)

myped = read.table("ped.txt", T)
Dped=myped[,-c(1:6)]
Dped=scale(Dped,T,T)

myphe=read.table("phe.txt", header=T)
model=glm(family=binomial, myphe[,3]~pc1)
