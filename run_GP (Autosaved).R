Names=c("pyrim","triazines","machine","housing","abalone")
N=c(100,300,500,700,1000)
nameind=1
Name=Names[nameind]
datname=sprintf("%s.data",Name)
domname=sprintf("%s.domain",Name)
Dat=read.table(file=datname,sep=",")
colnames(Dat)=read.table(file=domname,sep=":")$V1
source("./functions.R")
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]


m=N[nameind] #number of pairs in training data
tm=1000 #number of pairs in testing data

newDat=samplepairs(m,Dat,y)
tDat=samplepairs(tm,Dat,y)


norDat=newDat
ntDat0=tDat

cDat=convertDat(norDat)
#estimate fI
Kappa=4
Phi=0.5
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
		ctDat=convertDat(ntDat0)
		pre=pre_GP(getcov(ctDat$points,x,kernel1),Sig,fI)
		e1=loss_GP(pre,ctDat$pairs)/tm
		e2=loss_GP(fI,cDat$pairs)/m
		
		msg1=sprintf("kappa: %f, phi: %f",kap,phi)
		print(msg1)
		msg2=sprintf("testing err: %f, training err: %f",e1,e2)
		print(msg2)
		
		E1[l,j]=e1
		E2[l,j]=e2
		j=j+1
	}
	l=l+1
}
# savefile=sprintf("%s_GP%i.Rdata",Name,nameind)
# save(E1,E2,file=savefile)