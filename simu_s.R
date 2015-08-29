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
			SDat=spl(norDat,pow,k,inter=F,tDat=ntDat)
			sDat=SDat$sDat
			stDat=SDat$stDat
			result=SA_s(sDat,phi=2)
			b=result$b
			e1=loss_s(b,stDat)/tm
			e2=loss_s(b,sDat)/m
			msg=sprintf("p: %i, type: %i, m: %i ", p,t,m)
			print(msg)
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
save(TE1,TE2,file="./simu_s.Rdata")