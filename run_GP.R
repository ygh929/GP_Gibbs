Dat=read.table("pyrim.data",sep=",")
colnames(Dat)=read.table("pyrim.domain",sep=":")$V1
source("./functions.r")
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]


m=10 #number of pairs in training data
tm=20 #number of pairs in testing data

newDat=samplepairs(m,Dat,y)
tDat=samplepairs(tm,Dat,y)

#first norm to [-1,1]
nDat=normto1(newDat)
norDat=nDat$normDat
ranDat=nDat$ran
ntDat0=normbyran(tDat,ranDat)

cDat=convertDat(norDat)
#estimate fI
Kappa=1:10
Phi=(1:10)*0.2
E1=matrix(NA,length(Kappa),length(Phi))
E2=E1
E3=E1
l=1
for (kap in Kappa){
	kernel1<-function(x1,x2,kappa=kap){
		#kernal function to generate Sigma
		k=sum((x1-x2)^2)
		K=exp(-kappa*k/2)
		K
	}
	j=1
	for (phi in Phi){
		result=SA_GP(cDat,phi)
		fI=result$fI
		Sig=result$Sig
		#predict for test set and calculate error rate
		x=cDat$points
		e1=geterr_GP(ntDat0,x,Sig,fI)
		e2=geterr_GP(norDat,x,Sig,fI)
		e3=loss_GP(fI,cDat$pairs)/m
		
		msg1=sprintf("kappa: %f, phi: %f",kap,phi)
		print(msg1)
		msg2=sprintf("testing err: %f, training err: %f, loss at fI: %f",e1,e2,e3)
		print(msg2)
		
		E1[l,j]=e1
		E2[l,j]=e2
		E3[l,j]=e3
		j=j+1
	}
	l=l+1
}
save(E1,E2,E3,file="errors_GP.Rdata")