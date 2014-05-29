Dat=read.table("pyrim.data",sep=",")
colnames(Dat)=read.table("pyrim.domain",sep=":")$V1
source("./functions.r")
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]


m=100 #number of pairs in training data
tm=200 #number of pairs in testing data

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
for (kap in Kappa){
	kernel1<-function(x1,x2,kappa=kap){
		#kernal function to generate Sigma
		k=sum((x1-x2)^2)
		K=exp(-kappa*k/2)
		K
	}
	for (phi in Phi){
		result=SA_GP(cDat,phi)
		fI=result$fI
		Sig=result$Sig
		#predict for test set and calculate error rate
		err=rep(NA,tm)
		x=cDat$points
		#normalize test data first
		ntDat=ntDat0		
		for (i in 1:tm){
			xt=rbind(ntDat$X1[i,],ntDat$X2[i,])
			coSig=getcov(xt,x,kernel1)
			Y=pre_GP(coSig,Sig,fI)
			err[i]=I(Y[1]<Y[2])
		}
		e1=mean(err)
		err=rep(NA,m)
		ntDat=norDat
		for (i in 1:m){
			xt=rbind(ntDat$X1[i,],ntDat$X2[i,])
			coSig=getcov(xt,x,kernel1)
			Y=pre_GP(coSig,Sig,fI)
			err[i]=I(Y[1]<Y[2])
		}
		e2=mean(err)
		e3=loss_GP(fI,cDat$pairs)/m
		msg1=sprintf("kappa: %f, phi: %f",kap,phi)
		print(msg1)
		msg2=sprintf("testing err: %f, training err: %f, loss at fI: %f",e1,e2,e3)
		print(msg2)
		
				
	}
}
