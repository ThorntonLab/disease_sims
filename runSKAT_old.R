library(SKAT)

a=commandArgs(trailing=TRUE)

ccfile=a[1]
ofile=a[2]

raw=read.table(ccfile)

MBweights=function(x,gwas=FALSE,mfreq=0.05,ncontrols=3000)
{
	q = (sum(x[1:ncontrols])+1)/(2*ncontrols + 2)
	w=1/sqrt(ncontrols*q*(1-q))
	if( gwas == TRUE & q < mfreq )
	{
		w=0
	}
	return(w)
}

onchip = function(x,minfreq=0.05)
{
	q = sum(x)/length(x)
	rv=1
	if(q < minfreq)
	{
	rv=0
	}
	return(rv)
}

fet.p = function(x,ncontrols=3000)
{
	minor.controls = sum(x[1:ncontrols])
	major.controls = 2*ncontrols - minor.controls
	ncases=length(x)-ncontrols
	minor.cases = sum(x[(ncontrols+1):length(x)])
	major.cases = 2*ncases - minor.cases
	return ( fisher.test( matrix(c(minor.controls,major.controls,minor.cases,major.cases),nr=2) )$p.value )
}

data=as.matrix( raw[,1:(ncol(raw)-4)] + 1)

rm(raw)

data.mbweights=apply(data,2,MBweights)
data.mbweights.gwas=apply(data,2,MBweights,TRUE)
data.onchip = apply(data,2,onchip)

#3000 cases, followed by 3000 controls
b=c(rep(0,3000),rep(1,3000))
obj = SKAT_Null_Model( b ~ 1 )

MB.optimal.p = as.numeric(SKAT(data,obj,method="optimal",kernel="linear.weighted",weights=data.mbweights)$p.value)
MB.gwas.optimal.p = as.numeric(SKAT(data,obj,method="optimal",kernel="linear.weighted",weights=data.mbweights.gwas)$p.value)
SKAT.p = as.numeric(SKAT(data,obj)$p.value)
SKAT.gwas.p = as.numeric(SKAT(data[,which(data.onchip==1)],obj)$p.value)
SKAT.optimal.p = as.numeric(SKAT(data,obj,method="optimal")$p.value)
SKAT.gwas.optimal.p = as.numeric(SKAT(data[,which(data.onchip==1)],obj,method="optimal")$p.value)
MB.p = as.numeric(SKAT(data,obj,kernel="linear.weighted",weights=data.mbweights)$p.value)
MB.gwas.p = as.numeric(SKAT(data,obj,kernel="linear.weighted",weights=data.mbweights.gwas)$p.value)



#now, do the analog to our "esm" test
FET = apply(data,2,fet.p)
FET.rank = rank(FET,ties="random")
FET.gwas = FET
FET.gwas[which(data.onchip==0)]=50 #turn them into an impossibly-high p-value
FET.gwas.rank = rank(FET.gwas,ties="random")

#look @ 25 "best" markers
K=25
esm.p = as.numeric( SKAT(data[,which(FET.rank <= K)],obj)$p.value )
esm.gwas.p = as.numeric( SKAT(data[,which(FET.gwas.rank <= K)],obj)$p.value )
esm.optimal.p = as.numeric( SKAT(data[,which(FET.rank <= K)],obj,method="optimal")$p.value )
esm.gwas.optimal.p = as.numeric( SKAT(data[,which(FET.gwas.rank <= K)],obj,method="optimal")$p.value )

write.table(cbind(SKAT.p,SKAT.gwas.p,SKAT.optimal.p,SKAT.gwas.optimal.p,MB.p,MB.gwas.p,esm.p,esm.gwas.p,esm.optimal.p,esm.gwas.optimal.p,MB.optimal.p,MB.gwas.optimal.p),col.names=c("SKAT.p","SKAT.gwas.p","SKAT.optimal.p","SKAT.gwas.optimal.p","MB.p","MB.gwas.p","esm.p","esm.gwas.p","esm.optimal.p","esm.gwas.optimal.p","MB.optimal.p","MB.gwas.optimal.p"),file=ofile,row.names=FALSE,quote=FALSE)
