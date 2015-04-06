library(diseaseSims)
library(SKAT)

ccidx = read.table("ccindex.txt",colClasses=c("integer","numeric"))
ccdata = getCCblock("ccfile.bin.gz",procIndex(ccidx,2,0))
status = c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))


##For fun, let's only run SKAT on markers with MAF < 0.05 in controls
MAF = colSums(ccdata$genos[1:ccdata$ncontrols,])/(2*ccdata$ncontrols)

obj<-SKAT_Null_Model(status ~ 1, out_type="D")
skat.default=SKAT(ccdata$genos[,which(MAF<0.05)],obj)
skat.O=SKAT(ccdata$genos[,which(MAF<0.05)],obj,method="optimal.adj")
skat.linear=SKAT(ccdata$genos[,which(MAF<0.05)],obj,kernel="linear")
skat.linear.O=SKAT(ccdata$genos[,which(MAF<0.05)],obj,method="optimal.adj",kernel="linear")

print("SKAT p-values:")
print(paste("default: ",skat.default$p.value))
print(paste("O: ",skat.O$p.value))
print(paste("linear: ",skat.linear$p.value))
print(paste("linear/O: ",skat.linear.O$p.value))
