#sample pairs from the data

library(mvtnorm)#to generate rmvnorm

samplepairs<-function(m,Dat,y){
	#get sample of pairs from the data
	#values in X1 is preefered to values in X2
	nr=dim(Dat)[1]
	nc=dim(Dat)[2]
	newDat=list(X1=as.data.frame(NULL),X2=as.data.frame(NULL))

	fac.ind=sapply(Dat,is.factor)
	if (sum(fac.ind)>0){
		fDat=Dat[fac.ind]
		dummys=model.matrix(y~.,data=fDat)
		dummys=dummys[,2:ncol(dummys)]
		Dat=Dat[,-which(fac.ind==1)]
		Dat=cbind(Dat,dummys)
	}
	
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
	
SA_GP<-function(cDat,phi=0.5){
	#try to minimize q=l-log(p)
	kmax=1e5
	sig=0.5 #initial jump sd
	tolc=1 #initial tol
	x=cDat$points
	pairs=cDat$pairs
	Sig=getSigma(x,kernel1)
	n=nrow(x)
	r=sqrt(n)
	#initialize
	fI=rmvnorm(1,rep(0,n),Sig)
	fI=fI/sqrt(sum(fI^2))*r
	q=phi*loss_GP(fI,pairs)-dmvnorm(fI,mean=rep(0,n),sigma=Sig,log=TRUE)
	count=0
	#move
	for (k in 1:kmax){

		newfI=MoveonSph(fI,sig,r)
		newq=phi*loss_GP(newfI,pairs)-dmvnorm(newfI,mean=rep(0,n),sigma=Sig,log=TRUE)
	
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			fI=newfI
			q=newq
			count=count+1
		}
		if ((k%%1000)==0){
			#adjust tol and sig by accept rate
			accrate=count/1000
			if (accrate>0.2){
				tolc=tolc*5
			}
			if (accrate<0.2){
				sig=sig*0.75
			}
			if (accrate==0){
				break
			}
			count=0
		}
	}
	list(fI=fI,Sig=Sig,q=q)
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

geterr_GP<-function(ntDat,x,Sig,fI){
	tm=nrow(ntDat[[1]])
	err=rep(NA,tm)
	for (i in 1:tm){
			xt=rbind(ntDat$X1[i,],ntDat$X2[i,])
			coSig=getcov(xt,x,kernel1)
			Y=pre_GP(coSig,Sig,fI)
			err[i]=I(Y[1]<Y[2])
		}
	e=mean(err)
	e
}
normto1<-function(x){
	normDat=list(X1=NULL,X2=NULL)
	Tdata=rbind(x$X1,x$X2)
	ran=sapply(Tdata,range)
	for (j in 1:ncol(Tdata)){
		if (ran[1,j]!=ran[2,j]){
			Tdata[,j]=(Tdata[,j]-ran[1,j])/(ran[2,j]-ran[1,j])*2
		}else{
			Tdata[,j]=Tdata[,j]-ran[1,j]
		}
		
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
	for (j in 1:ncol(Tdat)){
		if (ran[1,j]!=ran[2,j]){
			Tdat[,j]=(Tdat[,j]-ran[1,j])/(ran[2,j]-ran[1,j])*2
		}else{
			Tdat[,j]=Tdat[,j]-ran[1,j]
		}
		
	}
	x$X1=Tdat[1:nr,]-1
	x$X2=Tdat[(nr+1):(2*nr),]-1
	x
}
spl<-function(Dat,c,k,tDat=NULL){
	#function to get cubic splines with k knots in each dimension
	x=rbind(Dat$X1,Dat$X2)
	m=nrow(Dat$X1)
	newx=as.data.frame(matrix(NA,nrow(x),0))
	if (!is.null(tDat)){
		m1=nrow(tDat$X1)
		xt=rbind(tDat$X1,tDat$X2)
		x=rbind(x,xt)
		newx=as.data.frame(matrix(NA,nrow(x),0))
	}
	
	for (j in 1:ncol(x)){
		tempx=x[,j]
		if (length(unique(tempx))<=2){
			names(tempx)=colnames(x)[j]
			newx=cbind(newx,tempx)
		}else{
			if (k>0){
				xi=min(tempx)+(range(tempx)[2]-range(tempx)[1])*(1:k)/(k+1)
				addx=matrix(NA,nrow(x),k+c)
				addx[,1:c]=outer(tempx,1:c,'^')
				h=outer(tempx,xi,'-')
				addx[,(c+1):(c+k)]=(h*(h>0))^c
				colnames(addx)=sprintf("%s%i",colnames(x)[j],1:(c+k))
				newx=cbind(newx,addx)	
			}else{
				addx=matrix(NA,nrow(x),k+c)
				addx[,1:c]=outer(tempx,1:c,'^')
				colnames(addx)=sprintf("%s%i",colnames(x)[j],1:(c+k))
				newx=cbind(newx,addx)	
				}
		}
	}
	rownames(newx)=rownames(x)
	sDat=list(X1=NULL,X2=NULL)
	sDat$X1=newx[1:m,]
	sDat$X2=newx[(m+1):(2*m),]
	output=list(sDat=sDat)
	if (!is.null(tDat)){
		stDat=list(X1=NULL,X2=NULL)
		stDat$X1=newx[(2*m+1):(2*m+m1),]
		stDat$X2=newx[(2*m+m1+1):(2*m+2*m1),]
		output=list(sDat=sDat,stDat=stDat)
	}
	output
}

loss_s=function(b,sDat){
	Y=matrix(NA,nrow=nrow(sDat[[1]]),ncol=2)
	for (i in 1:2){
		Y[,i]=as.matrix(sDat[[i]])%*%matrix(b,ncol=1)
	}
	loss=sum(Y[,1]<Y[,2])
	loss
}
	
SA_s<-function(sDat,phi=1){
	#try to minimize q=l-log(p)
	kmax=1e4
	sig=0.5 #initial jump sd
	tolc=1 #initial tol
	ps=ncol(sDat[[1]])
	r=sqrt(ps)
	Sig=diag(ps)
	#initialize
	b=rmvnorm(1,rep(0,ps),Sig)
	b=b/sqrt(sum(b^2))*r
	q=phi*loss_s(b,sDat)
	count=0
	#move
	for (k in 1:kmax){

		newb=MoveonSph(b,sig,r)
		newq=phi*loss_s(newb,sDat)
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			b=newb
			q=newq
			count=count+1
		}
		if ((k%%1000)==0){
			
			#adjust tol and sig by accept rate
			accrate=count/1000
			if (accrate>0.2){
				tolc=tolc*5
			}
			if (accrate<0.2){
				sig=sig*0.9
			}
			if (accrate==0){
				break
			}
			count=0
			print(k)
		}
	}
	list(b=b,q=q)
}


