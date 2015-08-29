source("./functions.R")
TE1=list(p3=NULL,p20=NULL)
TE2=list(p3=NULL,p20=NULL)
l=1
for (p in c(3,20)){
	E1=matrix(NA,5,2)
	E2=E1
	i=1
	for (t in 1:5){
		j=1
		for (m in c(100,500)){
			Gen=gendata(150,p,type=t)
			Dat=Gen$Dat
			nc=ncol(Dat)
			y=Dat[,nc]
			Dat=Dat[,1:(nc-1)]
			tm=2000 #number of pairs in testing data
			pow=3
			k=3
			newDat=samplepairs(m,Dat,y)
			tDat=samplepairs(tm,Dat,y)
			#sieve approach
			nDat=normto1(newDat)
			norDat=nDat$normDat
			ranDat=nDat$ran
			ntDat=normbyran(tDat,ranDat)	
			cDat=convertDat(norDat)	
			kap=2	
			phi=2
			kernel1=gen_ker(kap)
			result=SA_GP(cDat,phi)
			fI=result$fI
			Sig=result$Sig
			#predict for test set and calculate error rate
			x=cDat$points
			ctDat=convertDat(ntDat)
			pre=pre_GP(getcov(ctDat$points,x,kernel1),Sig,fI)
			e1=loss_GP(pre,ctDat$pairs)/tm
			e2=loss_GP(fI,cDat$pairs)/m
			msg=sprintf("p: %i, type: %i, m: %i ", p,t,m)
			print(msg)
			msg1=sprintf("kappa: %f, phi: %f",kap,phi)
			print(msg1)
			msg2=sprintf("testing err: %f, training err: %f",e1,e2)
			print(msg2)
			
			E1[i,j]=e1
			E2[i,j]=e2

			j=j+1

		}
		i=i+1
	}
	TE1[[l]]=E1
	TE2[[l]]=E2	
	l=l+1
}
save(TE1,TE2,file=sprintf("./simu_GPn_kp%i.Rdata",kap))