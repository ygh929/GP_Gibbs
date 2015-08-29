#sample pairs from the data

library(mvtnorm)#to generate rmvnorm
library(boot) #to use inv.logit
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

kernel1<-function(x1,x2,kappa=10){
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
	
SA_GP<-function(cDat,phi=0.5,fI0=NULL,kmax=5e4){
	
	#try to minimize q=l-log(p)
	sig=0.1 #initial jump sd
	tolc=1 #initial tol
	x=cDat$points
	pairs=cDat$pairs
	Sig=getSigma(x,kernel1)
	n=nrow(x)
	r=5
	#initialize
	if (is.null(fI0)){
		fI=rmvnorm(1,rep(0,n),Sig)
		fI=fI/sqrt(sum(fI^2))*r
	}else{
		fI=fI0
	}
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
			if (accrate>0.1){
				tolc=tolc*5
			}
			if (accrate<0.1){
				sig=sig*0.8
			}
			if (accrate==0){
				#break
			}
			count=0
		}
	}
	list(fI=fI,Sig=Sig,q=q)
}

lp_kappa<-function(kappa){
	d=(-2-1)*log(kappa)-2/kappa
	d	
}

gen_ker<-function(kappa=1){
	kernel1<-function(x1,x2){
		#kernal function to generate Sigma
		k=sum((x1-x2)^2)
		K=exp(-kappa*k/2)
		K
	}
	kernel1
}

SA_GP_k<-function(cDat,phi=1,fI0=NULL){
	#try to minimize q=l-log(p)
	kmax=1e5
	sig=list(b=0.1,kappa=0.1) #initial jump sd
	tolc=1 #initial tol
	x=cDat$points
	pairs=cDat$pairs
	Sig=getSigma(x,kernel1)
	n=nrow(x)
	r=5
	#initialize
	kappa=1
	if (is.null(fI0)){
		fI=rmvnorm(1,rep(0,n),Sig)
		fI=fI/sqrt(sum(fI^2))*r
	}else{
		fI=fI0
	}
	q=phi*loss_GP(fI,pairs)-dmvnorm(fI,sigma=Sig,log=TRUE)-lp_kappa(kappa)
	count=list(b=0,kappa=0)
	#move
	for (k in 1:kmax){
		newfI=MoveonSph(fI,sig$b,r)
		newq=phi*loss_GP(newfI,pairs)-dmvnorm(newfI,sigma=Sig,log=TRUE)-lp_kappa(kappa)
	
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			fI=newfI
			q=newq
			count$b=count$b+1
		}
		
		newkappa=abs(kappa+rnorm(1,0,sd=sig$kappa))
		newkernel=gen_ker(newkappa)
		newSig=getSigma(x,newkernel)
		newq=phi*loss_GP(fI,pairs)-dmvnorm(fI,sigma=newSig,log=TRUE)-lp_kappa(newkappa)
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			kappa=newkappa
			Sig=newSig
			q=newq
			count$kappa=count$kappa+1
		}
		if ((k%%1000)==0){
			
			
			#adjust tol and sig by accept rate
			accrate=list(b=count$b/1000,kappa=count$kappa/1000)
			for (i in 1:2){
				if (accrate[[i]]>0.1){
					tolc=tolc*2
				}
				if (accrate[[i]]<0.1){
					sig[[i]]=sig[[i]]*0.8
				}
				count[[i]]=0				
			}
		}
	}
	list(fI=fI,Sig=Sig,q=q,kappa=kappa)
}



getcov<-function(xt,x,ker){
	#x is all the observed n (we use 2m here) locations and xt is the new pair
	n=nrow(x)
	nt=nrow(xt)
	Sigma=matrix(1,n,nt)
	for (i in 1:n){
		for (j in 1:nt){
			Sigma[i,j]=ker(x[i,],xt[j,])
		}
	}
	Sigma
}


pre_GP<-function(coSig,Sig,f){
	
	Y=t(coSig)%*%solve(Sig)%*%matrix(f,ncol=1)
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
spl<-function(Dat,c,k,inter=FALSE,tDat=NULL){
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
			newx=cbind(newx,tempx)
			colnames(newx)[ncol(newx)]=colnames(x)[j]
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
		if (inter==TRUE && j!=ncol(x)){
			addx1=matrix(NA,nrow(x),0)
			for (l in (j+1):ncol(x)){
				tempadd=x[,j]*x[,l]
				if (length(unique(tempadd))>=2){
					addx1=cbind(addx1,tempadd)
					colnames(addx1)[ncol(addx1)]=sprintf("%s.%s",colnames(x)[j],colnames(x)[l])
					#colnames(addx1)=sprintf("%s.%s",colnames(x[,j]),colnames(x[,l]))
				}
				
			}
			newx=cbind(newx,addx1)
			
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
	
SA_s<-function(sDat,phi=1,b0=NULL){
	#try to minimize q=l-log(p)
	kmax=1e5
	sig=0.5 #initial jump sd
	tolc=1 #initial tol
	ps=ncol(sDat[[1]])
	r=sqrt(ps)
	Sig=diag(ps)
	#initialize
	if (is.null(b0)){
		b=rmvnorm(1,rep(0,ps),Sig)
		b=b/sqrt(sum(b^2))*r
	}else{
		b=b0
	}
	
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
				sig=sig*0.5
			}
			if (accrate==0){
				#break
			}
			count=0
			
		}
	}
	list(b=b,q=q)
}

I1<-function(z,sig=0.01){
	Ind=inv.logit(z/sig)
	Ind
}

loss_s1=function(b,sDat){
	Y=matrix(NA,nrow=nrow(sDat[[1]]),ncol=2)
	for (i in 1:2){
		Y[,i]=as.matrix(sDat[[i]])%*%matrix(b,ncol=1)
	}
	loss=sum(I1(Y[,2]-Y[,1]))
	loss
}


SA_s1<-function(sDat,phi=1){
	#try to minimize q=l-log(p)
	kmax=5e4
	sig=0.5 #initial jump sd
	tolc=1 #initial tol
	ps=ncol(sDat[[1]])
	Sig=diag(ps)
	#initialize
	b=rmvnorm(1,rep(0,ps),Sig)
	q=phi*loss_s1(b,sDat)-dmvnorm(b,log=T)
	count=0
	#move
	for (k in 1:kmax){

		newb=b+rnorm(ps,0,sig)
		newq=phi*loss_s1(newb,sDat)-dmvnorm(newb,log=T)
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			b=newb
			q=newq
			count=count+1
		}
		if ((k%%1000)==0){

			#adjust tol and sig by accept rate
			accrate=count/1000
			msg=sprintf("k=%i, accrate=%f, current err=%i",k,accrate,loss_s(b,sDat))
			print(msg)
			if (accrate>0.2){
				tolc=tolc*5
			}
			if (accrate<0.2){
				sig=sig*0.9
			}
			if (accrate==0){
				#break
			}
			count=0
		}
	}
	list(b=b,q=q)
}


SA_vs<-function(sDat,phi=1,rho=0.5,b0=NULL){
	#try to minimize q=l-log(p)
	kmax=1e5
	sig=0.5 #initial jump sd
	tolc=1 #initial tol
	ps=ncol(sDat[[1]])
	r=sqrt(ps*rho)
	#initialize
	if (is.null(b0)){
		b=rmvnorm(1,rep(0,ps))
		S=rbinom(ps,1,rho)
		b=b*S
		b=b/sqrt(sum(b^2))*r
	}else{
		b=b0
		S=(b!=0)	
	}

	adj0=dnorm(0)+log(rho)-log(1-rho)
	bn0=b[S]
	nn0=length(bn0)
	q=phi*loss_s(b,sDat)+adj0*(ps-nn0)
	count=0
	#move
	for (k in 1:kmax){	
		move=0
		change=sample(1:ps,1)
		newS=S
		newS[change]=1-S[change]
		if (S[change]==1){
			newb=b*newS
			newb=newb/sqrt(sum(newb^2))*r
			newq=phi*loss_s(newb,sDat)+adj0*(ps-sum(newS))
		}else{
			newb=b
			newb[change]=rnorm(1,mean=0,sd=sig)
			newb=newb/sqrt(sum(newb^2))*r
			newq=phi*loss_s(newb,sDat)+adj0*(ps-sum(newS))
		}
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			b=newb
			S=newS
			q=newq
			move=move+1
		}	
		newb=b
		newb[S]=b[S]+rnorm(sum(S),0,sig)
		newb=newb/sqrt(sum(newb^2))*r
		newq=phi*loss_s(newb,sDat)+adj0*(ps-sum(newS))
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			b=newb
			q=newq
			move=move+1
		}
		if (move>0){
			count=count+1
		}	
		if ((k%%1000)==0){

			#adjust tol and sig by accept rate
			accrate=count/1000
			msg=sprintf("k=%i, accrate=%f, current err=%i",k,accrate,loss_s(b,sDat))
			print(msg)
			if (accrate>0.2){
				tolc=tolc*5
			}
			if (accrate<0.2){
				sig=sig*0.95
			}
			if (accrate==0){
				#break
			}
			count=0
			
		}
	}
	list(b=b,q=q)
}

gendata<-function(n,p,type=1,errsd=0){
	b=rnorm(p)
	X=matrix(runif(n*p)*10,n,p)
	y1<-function(x){
		x%*%b
	}
	y2<-function(x){
		x^4%*%b
	}
	y3<-function(x){
		exp(x)%*%b
	}
	y4<-function(x){
		(x%*%b)^2
	}
	y5<-function(x){
		t(rmvnorm(1,sigma=getSigma(x,kernel1)))
	}
	y=list(y1,y2,y3,y4,y5)
	Y=y[[type]](X)+rnorm(n,0,errsd)
	Dat=cbind(X,Y)
	colnames(Dat)=1:ncol(Dat)
	list(Dat=Dat,b=b)
}
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

SA_BLP<-function(cDat,vcDat){
	vZ=GetXZ(vcDat)$Z
	kmax=3e3
	tolc=1
	sig_kappa=0.2
	sig_sig=0.2
	count=0
	for (k in 1:kmax){
		if (k==1){
			#initialize
			kappa=2
			sig=1
			ker=gen_ker(kappa)
			Sig=getSigma(cDat$points,ker)
			ParC=PartialCov(cDat$points,vcDat$points,ker)
			q=sum(pre_BLP(vcDat$pairs,Sig,ParC,X,Z,sig=sig)!=vZ)

		}
		newkappa=abs(kappa+rnorm(1,0,sig_kappa))
		newsig=abs(sig+rnorm(1,0,sig_sig))
		newker=gen_ker(newkappa)
		newSig=getSigma(cDat$points,newker)
		newParC=PartialCov(cDat$points,vcDat$points,newker)
		possibleErr=tryCatch(
			{
				newq=sum(pre_BLP(vcDat$pairs,newSig,newParC,X,Z,sig=newsig)!=vZ)
				newq
			},
			error=function(e)e
		)
		if (inherits(possibleErr,"error")){
			next
		}
		
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			sig=newsig
			kappa=newkappa
			q=newq
			count=count+1
		}
		if ((k%%200)==0){

			#adjust tol and sig by accept rate
			accrate=count/200

			if (accrate>0.2){
				tolc=tolc*5
			}
			if (accrate<0.2){
				sig_sig=sig_sig*0.8
				sig_kappa=sig_kappa*0.8
			}
			if (accrate==0){
				break
			}
			count=0
		}
	}
	list(kappa=kappa,sig=sig,err=q/length(vZ))
}

SA_GP_v<-function(vcDat,fI0,cDat,kappa0=10,kmax=2e3){
	sig=list(kappa=5) #initial jump sd
	tolc=5 #initial tol
	x=vcDat$points
	x0=cDat$points
	pairs=vcDat$pairs
	kappa=kappa0
	ker=gen_ker(kappa)
	Sig=getSigma(x0,ker)
	coSig=getcov(x,x0,ker)
	n=nrow(x)
	#initialize

	fI=fI0
	q=loss_GP(pre_GP(coSig,Sig,fI),pairs)
	count=list(kappa=0)
	#move
	for (k in 1:kmax){
		newkappa=abs(kappa+rnorm(1,0,sd=sig$kappa))
		newker=gen_ker(newkappa)
		newSig=getSigma(x0,newker)
		newcoSig=getcov(x,x0,newker)
		possibleErr=tryCatch(
		{newq=loss_GP(pre_GP(newcoSig,newSig,fI),pairs)
			newq},
		error=function(e)e
		)
		if (inherits(possibleErr,"Error")){
			next
		}
		temp=runif(1,0,1)
		if (temp<(exp(tolc*(q-newq)))){
			kappa=newkappa
			q=newq
			count$kappa=count$kappa+1
			Sig=newSig
		}
		if ((k%%100)==0){
			
			
			#adjust tol and sig by accept rate
			accrate=list(kappa=count$kappa/100)
		
			if (accrate$kappa>0.2){
				tolc=tolc*5
			}
			if (accrate$kappa<=0.2){
				sig$kappa=sig$kappa*0.5
			}
			if (accrate$kappa==0){
				break
			}
			count$kappa=0
						
			
		}
	}
	list(kappa=kappa,q,fI=fI,Sig=Sig)
}