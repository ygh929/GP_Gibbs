#sample pairs from the data
Dat=read.table("abalone.data",sep=",")
colnames(Dat)=read.table("abalone.domain",sep=":")$V1

samplepairs<-function(m,Dat,y){
	nr=dim(Dat)[1]
	nc=dim(Dat)[2]
	newDat=list(X1=as.data.frame(NULL),X2=as.data.frame(NULL))
	colnames(newDat$X1)=colnames(Dat)
	colnames(newDat$X2)=colnames(Dat)
	fac.ind=sapply(Dat,is.factor)
	fDat=Dat[fac.ind]
	dummys=model.matrix(y~.,data=fDat)
	dummys=dummys[,2:ncol(dummys)]
	Dat=Dat[,-which(fac.ind==1)]
	Dat=cbind(Dat,dummys)
	i=1
	while (i<=m){
		i1=sample(1:nr,1)
		i2=sample(1:nr,1)
		if (y[i1]>y[i2]){
			newDat$X1=rbind(newDat$X1,Dat[i1,])
			newDat$X2=rbind(newDat$X2,Dat[i2,])
			i=i+1
		}
		if (y[i1]<y[i2]){
			newDat$X1=rbind(newDat$X1,Dat[i2,])
			newDat$X2=rbind(newDat$X2,Dat[i1,])
			i=i+1			
		}
	}
	newDat
}

kernel1<-function(x1,x2,kappa=1){
	k=sum((x1-x2)^2)
	K=exp(-kappa*k/2)
	K
}

cubicspl<-function(x,k){
	newx=as.data.frame(matrix(NA,nrow(x),0))
	for (j in 1:ncol(x)){
		tempx=x[,j]
		if (length(unique(tempx))<=2){
			names(tempx)=colnames(x)[j]
			newx=cbind(newx,tempx)
			
		}else{
			xi=min(tempx)+(range(tempx)[2]-range(tempx)[1])*(1:k)/(k+1)
			addx=matrix(NA,nrow(x),k+3)
			addx[,1:3]=outer(tempx,1:3,'^')
			h=outer(tempx,xi,'-')
			addx[,4:(3+k)]=(h*(h>0))^3
			colnames(addx)=sprintf("%s%i",colnames(x)[j],1:(k+3))
			newx=cbind(newx,addx)	
		}
	}
	newx
}

