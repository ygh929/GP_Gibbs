GetXZ<-function(cDat){
	points=cDat$points
	pairs=cDat$pairs
	m=dim(pairs)[1]
	n=dim(points)[1]
	X=matrix(0,nrow=m,ncol=n)
	Z=matrix(NA,nrow=m,ncol=1)
	for (i in 1:m){
		X[i,min(pairs[i])]=1
		X[i,max(pairs[i])]=-1
		Z[i,1]=2*I(pairs[i,1]>pairs[i,2])-1
	}
	list(X,Z)
}
#set the size of train and validate data set
n1=50
n2=50
#generate training data and calculate X,Z
Gen1=gendata(n1,p,type=1)
Dat=Gen1$Dat
cDat=convertDat(Dat)
M=getXZ(cDat)
X=M$X
Z=M$Z
Sig=getSigma()
