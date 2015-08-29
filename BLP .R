source("./functions.R")
GetXZ<-function(cDat){
	points=cDat$points
	pairs=cDat$pairs
	m=dim(pairs)[1]
	n=dim(points)[1]
	X=matrix(0,nrow=m,ncol=n)
	Z=matrix(NA,nrow=m,ncol=1)
	for (i in 1:m){
		X[i,min(pairs[i,])]=1
		X[i,max(pairs[i,])]=-1
		Z[i,1]=2*I(pairs[i,1]<pairs[i,2])-1
	}
	list(X=X,Z=Z)
}

PartialCov<-function(X0,X1,ker){
	m=nrow(X0)
	n=nrow(X1)
	ParCov=matrix(1,nrow=n,ncol=m)
	for (i in 1:m){
		for (j in 1:n){
			ParCov[j,i]=ker(X0[i,],X1[j,])
		}
	}
	ParCov
}

pre_BLP<-function(tpairs,Sig,ParC,X,Z,sig=1){
	n1=nrow(Sig)
	n2=nrow(tpairs)
	PC=matrix(NA,n2,n1)
	for (i in 1:n2){
		PC[i,]=ParC[tpairs[i,1],]-ParC[tpairs[i,2],]
	}
	invSig=solve(Sig)
	out=sign(PC%*%invSig%*%solve(sig^2*invSig+t(X)%*%X)%*%t(X)%*%Z)
	out
}
m1=25
m2=75
nt=200
Gen=gendata(150,5,type=1)
Dat=Gen$Dat
nc=ncol(Dat)
y=Dat[,nc]
Dat=Dat[,1:(nc-1)]
#generate training data and calculate X,Z
newDat=samplepairs(m1,Dat,y)
valDat=samplepairs(m2,Dat,y)
tesDat=samplepairs(nt,Dat,y)
cDat=convertDat(newDat)
vcDat=convertDat(valDat)
tcDat=convertDat(tesDat)
tZ=GetXZ(tcDat)$Z
n1=nrow(cDat$points)
n2=nrow(vcDat$points)
M=GetXZ(cDat)
X=M$X
Z=M$Z

result=SA_BLP(cDat,vcDat)
ker=gen_ker(result$kappa)
Sig=getSigma(cDat$points,ker)
parC=PartialCov(cDat$points,tcDat$points,ker)
err=sum(pre_BLP(tcDat$pairs,Sig,parC,X,Z,sig=result$sig)!=tZ)/nt
print(err)
