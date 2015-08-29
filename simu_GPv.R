source("./functions.R")
TE1=list(p3=NULL,p20=NULL)
TE2=list(p3=NULL,p20=NULL)
TE3=list(p3=NULL,p20=NULL)
mt=200
l=1
for (p in c(3,20)){
	E1=matrix(NA,5,2)
	E2=E1
	E3=E1
	i=1
	for (t in 1:5){
		j=1
		for (m in c(100,500)){
			m1=1
			m2=m-m1
			Gen=gendata(150,p,type=t)
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
			result0=SA_GP(cDat,2,kmax=1e4)
			result=SA_GP_v(vcDat,result0$fI,cDat,kmax=1e3)
			x=cDat$points
			kappa=result$kappa
			ker=gen_ker(kappa)
			fI=result$fI
			Sig=result$Sig
			pre=pre_GP(getcov(tcDat$points,x,ker),Sig,fI)
			e1=loss_GP(pre,tcDat$pairs)/mt
			e2=loss_GP(pre_GP(getcov(vcDat$points,x,ker),Sig,fI),vcDat$pairs)/m2
			e3=loss_GP(fI,cDat$pairs)/m1
			msg=sprintf("p: %i, type: %i, m: %i ", p,t,m)
			print(msg)
			msg1=sprintf("kappa: %f",result$kappa)
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

save(TE1,TE2,TE3,file="./simu_GPv.Rdata")
