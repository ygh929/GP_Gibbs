Dat=read.table("pyrim.data",sep=",")
colnames(Dat)=read.table("pyrim.domain",sep=":")$V1
source("./functions.r")
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]


m=100 #number of pairs in training data
tm=100 #number of pairs in testing data
N=100# number of samples in posterior sampling

newDat=samplepairs(m,Dat,y)
tDat=samplepairs(tm,Dat,y)

#first norm to [-1,1]
nDat=normto1(newDat)
norDat=nDat$normDat
ranDat=nDat$ran
cDat=convertDat(norDat)
#estimate fI
result=SA_GP(cDat)
fI=result$fI
Sig=result$Sig
#predict for test set and calculate error rate
err=rep(NA,tm)
x=cDat$points
#normalize test data first
ntDat=normbyran(tDat,ranDat)
for (i in 1:tm){
	xt=rbind(ntDat$X1[i,],ntDat$X2[i,])
	coSig=getcov(xt,x,kernel1)
	Y=pre_GP(coSig,Sig,fI)
	err[i]=I(Y[1]<Y[2])
}
mean(err)