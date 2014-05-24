Dat=read.table("abalone.data",sep=",")
colnames(Dat)=read.table("abalone.domain",sep=":")$V1
source("./functions.r")
y=Dat[,9]
Dat=Dat[,1:8]


m=10 #number of pairs in training data
tm=100 #number of pairs in testing data
N=100# number of samples in posterior sampling

newDat=samplepairs(m,Dat,y)
tDat=samplepairs(tm,Dat,y)
#GP approach_sampling
x=rbind(newDat$X1,newDat$X2)
Sig=getSigma(x,kernel1)
p=ncol(newDat$X1)
post=list(Y1=NULL,Y2=NULL)
est=list(Y1=1:m,Y2=1:m)
i=1
while (i<=N){
	Y=mvrnorm(1,rep(0,2*m),Sig)
	est$Y1=Y[1:m]
	est$Y2=Y[(m+1):(2*m)]
	S=losssum(est)
	temp=runif(1,0,1)
	if (temp<(exp(-S))){
		post$Y1=cbind(post$Y1,est$Y1)
		post$Y2=cbind(post$Y2,est$Y2)
		i=i+1
	}
}
post_GP=post

err=(1:tm)*0
for (i in 1:tm){
	xt=rbind(tDat$X1[i,],tDat$X2[i,])
	coSig=getcov(xt,x,kernel1)
	tempY=matrix(NA,nrow=2,ncol=ncol(post$Y1))
	
	for (j in 1:ncol(post$Y1)){
		f=cbind(post$Y1[,j],post$Y2[,j])
		tempY[,j]=pre_GP(coSig,Sig,f)
	}
	err[i]=I(mean(tempY[1,])<mean(tempY[2,]))
	
}


#sieve approach
nto1=normto1(newDat)
normDat=nto1$normDat
ran=nto1$normDat
sieveDat=list(X1=NULL,X2=NULL)
sieveDat$X1=cubicspl(normDat$X1,1)
sieveDat$X2=cubicspl(normDat$X2,1)
ps=ncol(sieveDat$X1)
postb=NULL
est=list(Y1=1:m,Y2=1:m)
i=1
while(i<=N){
	b=matrix(rnorm(ps,0,1),nrow=ps,ncol=1)
	est$Y1=as.matrix(sieveDat$X1) %*% b
	est$Y2=as.matrix(sieveDat$X2) %*% b
		S=losssum(est)
	temp=runif(1,0,1)
	if (temp<(exp(-S))){
		postb=cbind(postb,b)
		i=i+1
	}
}

