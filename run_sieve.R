Dat=read.table("pyrim.data",sep=",")
colnames(Dat)=read.table("pyrim.domain",sep=":")$V1
source("./functions.r")
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]

m=100 #number of pairs in training data
tm=2000 #number of pairs in testing data
newDat=samplepairs(m,Dat,y)
tDat=samplepairs(tm,Dat,y)
#sieve approach
nDat=normto1(newDat)
norDat=nDat$normDat
ranDat=nDat$ran
ntDat=normbyran(tDat,ranDat)

Pow=1:3
K=0:2
i=1
j=1
E1=matrix(NA,nrow=length(Pow),ncol=length(K))
E2=E1

for (pow in Pow){
	j=1
	for (k in K){
		SDat=spl(norDat,pow,k,ntDat)
		sDat=SDat$sDat
		stDat=SDat$stDat
		result=SA_s(sDat,1)
		b=result$b
		e1=loss_s(b,stDat)/tm
		e2=loss_s(b,sDat)/m
		msg1=sprintf("power: %f, nknots: %f",pow,k)
		print(msg1)
		msg2=sprintf("testing err: %f, training err: %f",e1,e2)
		print(msg2)
		E1[i,j]=e1
		E2[i,j]=e2
		j=j+1
	}
	i=i+1
}
save(E1,E2,file="errors_s.Rdata")



