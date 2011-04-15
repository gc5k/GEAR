a=2
library("batch")
parseCommandArgs()

for(idx in 0:a) {
  pedfile = paste(idx, "ped.txt", sep="")
  phefile = paste(idx, "phe.txt", sep="")
  scorefile = paste(idx, "score.txt", sep="")
  mydata=read.table(pedfile, header=T)
  UIdx=which(mydata[,3]==0)
  RIdx=which(mydata[,3]!=0)
  d=mydata[,-c(1:6)]
  D=scale(d, T, T)

  UD=D[UIdx,]
  pca=prcomp(UD)
  plot(scale(pca$x[,1],T,T), scale(pca$x[,2],T,T))

  kc1=D%*%pca$rotation[,1]
  kc2=D%*%pca$rotation[,2]
  plot(scale(kc1,T,T), scale(kc2,T,T))

  myphe=as.matrix(read.table(phefile, header=T))
  model=glm(family=binomial, myphe[,3]~kc1+kc2)
  plot(myphe[,3]-model$fitted.values)
  
  scr = matrix(0, nrow(myphe), 4)
  scr[,1:3]=myphe[,1:3]
  scr[,4]=myphe[,3]-model$fitted.values
  
  title="FID PID Status Score";
  out=file(scorefile, "w")
  write(title, file=out, ncolumns=ncol(scr))
  for(i in 1:nrow(scr)) {
    write(scr[i,], file=out, ncolumns=ncol(scr))
  }
  close(out)
}
