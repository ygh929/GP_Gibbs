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
#sieve approach
nDat=normto1(newDat)
norDat=nDat$normDat
ranDat=nDat$normDat
ntDat=normbyran(tDat,ranDat)

sDat=cubicspl(norDat,1)
stDat=cubicspl(ntDat,1)
b=SA_s(sDat,1)
loss_s(b,stDat)

