Names=c("pyrim","triazines","machine","housing","abalone")
N=c(100,300,500,700,1000)
nameind=2
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
tm=2000 #number of pairs in testing data
Pow=3:3
K=2:2

TE1=list(NA)
TE2=list(NA)
E1=matrix(NA,nrow=length(Pow),ncol=length(K))
E2=E1
for (ite in 1:1){
	newDat=samplepairs(m,Dat,y)
	tDat=samplepairs(tm,Dat,y)
	#sieve approach
	nDat=normto1(newDat)
	norDat=nDat$normDat
	ranDat=nDat$ran
	ntDat=normbyran(tDat,ranDat)
	
	i=1
	for (pow in Pow){
		j=1
		for (k in K){
			SDat=spl(norDat,pow,k,inter=F,tDat=ntDat)
			sDat=SDat$sDat
			stDat=SDat$stDat
			result=SA_vs(sDat,phi=2,rho=0.4)
			b=result$b
			e1=loss_s(b,stDat)/tm
			e2=loss_s(b,sDat)/m
			msg=sprintf("it: %i", ite)
			print(msg)
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
	TE1[[ite]]=E1
	TE2[[ite]]=E2	
}

savefile=sprintf("%s_s%i.Rdata",Name,nameind)
save(TE1,TE2,file=savefile)



