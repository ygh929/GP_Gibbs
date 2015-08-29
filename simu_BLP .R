source("./functions.R")
TE1=list(p3=NULL,p20=NULL)
TE2=list(p3=NULL,p20=NULL)
TE3=list(p3=NULL,p20=NULL)
mt=2000
l=1
for (p in c(3,20)){
	E1=matrix(NA,5,2)
	E2=E1
	E3=E1
	i=1
	for (t in 1:5){
		j=1
		for (m in c(100,500)){
			m1=m/2
			m2=m-m1
			Gen=gendata(150,p,type=t,1)
			Dat=Gen$Dat
			nc=ncol(Dat)
			y=Dat[,nc]
			Dat=Dat[,1:(nc-1)]
			#generate training data and calculate X,Z
			newDat=samplepairs(m1,Dat,y)
			valDat=samplepairs(m2,Dat,y)
			tesDat=samplepairs(mt,Dat,y)
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
			e1=sum(pre_BLP(tcDat$pairs,Sig,parC,X,Z,sig=result$sig)!=tZ)/mt
			e2=result$err
			e3=sum(pre_BLP(cDat$pairs,Sig,Sig,X,Z,sig=result$sig)!=Z)/m1
			msg=sprintf("p: %i, type: %i, m: %i ", p,t,m)
			print(msg)
			msg1=sprintf("kappa: %f, phi: %f",result$kappa,result$sig)
			print(msg1)
			msg2=sprintf("testing err: %f, validating err: %f, training err: %f",e1,e2,e3)
			print(msg2)
			E1[i,j]=e1
			E2[i,j]=e2	
			E3[i,j]=e3
			j=j+1
		}
		i=i+1
	}
	TE1[[l]]=E1
	TE2[[l]]=E2
	TE3[[l]]=E3
	l=l+1
}

save(TE1,TE2,TE3,file="./simu_BLP.Rdata")
