library(diseaseSims)
library(buRden)

##Read index files
simidx = read.table("simindex.txt",colClasses=c("integer",rep("numeric",3)))
ccidx = read.table("ccindex.txt",colClasses=c("integer","numeric"))

##Read in the relefant data
esizes = getEsizes("effects.bin.gz",procIndex(simidx,2,0))
ccdata = getCCblock("ccfile.bin.gz",procIndex(ccidx,2,0))

##make a 0/1 list for controls/cases
status = c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))

##This is 2*(population size)
twoN=4e4

##Do logit regressions under various models
pvblock=makePVblock(ccdata,esizes,status,"additive")
pvblock[,"popfreq"]=pvblock[,"popfreq"]/twoN
##Add replicate id number (0) to data.frame:
pvblock = data.frame(cbind(repid=rep(0,nrow(pvblock)),pvblock))
##Request an access lock on the output file
OFILE="smarkerAdditive.txt.gz"
LOCK = Lock(OFILE)
OF=gzfile(OFILE,"a")

##Typically, you would do col.names=FALSE here, and ensure that the column
##names were written once and only once to the top of the file.
write.table(pvblock,OF,row.names=FALSE,append=TRUE,quote=FALSE,col.names=TRUE)
close(OF)

##Critical: failure to unlock will result in access lock persisting
##entire duration of R session
##alternate method will call the Rcpp destructors:
##rm(LOCK)
##gc()
Unlock(LOCK)

pvblock=makePVblock(ccdata,esizes,status,"recessive")
pvblock[,"popfreq"]=pvblock[,"popfreq"]/twoN
pvblock = data.frame(cbind(repid=rep(0,nrow(pvblock)),pvblock))

OFILE="smarkerRecessive.txt.gz"
LOCK = Lock(OFILE)
OF=gzfile(OFILE,"a")
write.table(pvblock,OF,row.names=FALSE,append=TRUE,quote=FALSE,col.names=TRUE)
close(OF)
Unlock(LOCK)

pvblock=makePVblock(ccdata,esizes,status,"dominant")
pvblock[,"popfreq"]=pvblock[,"popfreq"]/twoN
pvblock = data.frame(cbind(repid=rep(0,nrow(pvblock)),pvblock))
OFILE="smarkerDominant.txt.gz"
LOCK = Lock(OFILE)
OF=gzfile(OFILE,"a")
write.table(pvblock,OF,row.names=FALSE,append=TRUE,quote=FALSE,col.names=TRUE)
close(OF)
Unlock(LOCK)

