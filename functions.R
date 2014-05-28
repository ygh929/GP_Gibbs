#sample pairs from the data

library(mvtnorm)#to generate rmvnorm

samplepairs<-function(m,Dat,y){
	#get sample of pairs from the data
	#values in X1 is preefered to values in X2
	nr=dim(Dat)[1]
	nc=dim(Dat)[2]
	newDat=list(X1=as.data.frame(NULL),X2=as.data.frame(NULL))

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
	colnames(newDat$X1)=colnames(Dat)
	colnames(newDat$X2)=colnames(Dat)
	newDat
}

kernel1<-function(x1,x2,kappa=1){
	#kernal function to generate Sigma
	k=sum((x1-x2)^2)
	K=exp(-kappa*k/2)
	K
}

getSigma<-function(x,ker){
	#Get the covariance function given x values and kernel
	n=nrow(x)
	Sigma=matrix(1,n,n)
	for (i in 2:n){
		for (j in 1:(i-1)){
			Sigma[i,j]=ker(x[i,],x[j,])
			Sigma[j,i]=Sigma[i,j]
		}
	}	
	Sigma
}

MoveonSph=function(pos,sig,r){
	p=length(pos)
	newpos=pos+rnorm(p,0,sig)
	newpos=newpos/sqrt(sum(newpos^2))*r
	newpos
}

convertDat<-function(newDat){
	TDat=rbind(newDat$X1,newDat$X2)
	Points=unique(TDat)
	TDat=as.matrix(TDat)
	Points=as.matrix(Points)
	n=nrow(Points)
	pairs=matrix(NA,nrow=nrow(newDat$X1),ncol=2)
	for (i in 1:n){
		
		row.match=apply(TDat,1,identical,Points[i,])
		ind.match=which(row.match)
		pairs[ind.match]=i
	}
	list(points=Points,pairs=pairs)
}

loss_GP=function(fI,pairs){
	Y=pairs
	for (i in 1:length(Y)){
		Y[i]=fI[pairs[i]]
	}
	loss=sum(Y[,1]<Y[,2])
	loss
}
	
SA_GP<-function(cDat){
	#try to minimize q=l-log(p)
	kmax=1e4
	kappa=10
	sig=0.5 #initial jump sd
	tol=1e5 #initial tol
	x=cDat$points
	pairs=cDat$pairs
	Sig=getSigma(x,kernel1)
	n=nrow(x)
	r=sqrt(n)
	#initialize
	fI=rmvnorm(1,rep(0,n),Sig)
	fI=fI/sqrt(sum(fI^2))*r
	q=kappa*loss_GP(fI,pairs)#-dmvnorm(fI,mean=rep(0,n),sigma=Sig,log=TRUE)
	print(q)
	#move
	for (k in 1:kmax){
		count=0
		newfI=rmvnorm(1,rep(0,n),Sig)
		newFi=newfI/sqrt(sum(newfI^2))*r
		newq=kappa*loss_GP(newfI,pairs)#-dmvnorm(newfI,mean=rep(0,n),sigma=Sig,log=TRUE)
		print(newq)
		if ((newq-tol)<q){
			temp=runif(1,0,1)
			if (temp<(exp(q-newq))){
				fI=newfI
				q=newq
				count=count+1
			}
		}
		if (k%%100==0){
			#adjust tol and sig by accept rate
			accrate=count/100
			if (accrate>0.5){
				tol=tol/2
			}
			if (accrate<0.2){
				#sig=sig/2
			}
			if (accrate==0){
				#break
			}
			count=0
			print(accrate)
		}
	}
	list(fI=fI,Sig=Sig,L=loss_GP(fI,pairs))
}


getcov<-function(xt,x,ker){
	#x is all the observed n (we use 2m here) locations and xt is the new pair
	n=nrow(x)
	Sigma=matrix(1,n,2)
	for (i in 1:n){
		for (j in 1:2){
			Sigma[i,j]=ker(x[i,],xt[j,])
		}
	}
	Sigma
}




pre_GP<-function(coSig,Sig,f){
	
	Y=t(coSig)%*%Sig%*%matrix(f,ncol=1)
	Y
}
normto1<-function(x){
	normDat=list(X1=NULL,X2=NULL)
	Tdata=rbind(x$X1,x$X2)
	ran=sapply(Tdata,range)
	for (j in 1:ncol(Tdata)){
		Tdata[,j]=(Tdata[,j]-ran[1,j])/(ran[2,j]-ran[1,j])*2
	}
	m=nrow(x[[1]])
	normDat$X1=Tdata[1:m,]-1
	normDat$X2=Tdata[(m+1):(2*m),]-1
	list(normDat=normDat,ran=ran)
}
normbyran<-function(x,ran){
	Tdat=rbind(x$X1,x$X2)
	nr=nrow(x$X1)
	nc=ncol(Tdat)
	for (j in 1:nc){
		Tdat[,j]=(Tdat[,j]-ran[1,j])/(ran[2,j]-ran[1,j])*2-1
	}
	x$X1=Tdat[1:nr,]
	x$X2=Tdat[(nr+1):(2*nr),]
	x
}
cubicspl<-function(x,k){
	#function to get cubic splines with k knots in each dimension
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
	rownames(newx)=rownames(x)
	newx
}

losssum<-function(est){
	#est has two elements in the list $Y1 and $Y2
	S=sum(est$Y2>est$Y1)
	S
}


