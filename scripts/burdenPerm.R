library(diseaseSims)
library(buRden)

n=commandArgs(trailing=TRUE)
## Replicate ID
repID = as.integer(n[1])
CCIDX=n[2]
CCFILE=n[3]
## output file name
OFILENAME = n[4]
## How many perms to do
NPERMS = as.integer(n[5])

ccidx = read.table(CCIDX,colClasses=c("integer","numeric"))
ccdata = getCCblock(CCFILE,procIndex(ccidx,2,repID))
status = c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))

##Minor allele frequencies in controls
MAF = colSums(ccdata$genos[1:ccdata$ncontrols,])/(2*ccdata$ncontrols)

##Now, permute all burden tests NPERMS times

##Additional params

##number of markers for the ESM test
K=50

LLCMAF = 0.05

perms = allBurdenStatsPerm(ccdata$genos[,which(MAF < 0.05)],status,NPERMS,K,LLCMAF)

##Obtain a POSIX file lock on the output file
LOCK = Lock(OFILENAME)

F = gzfile(OFILENAME,"a")
##The output order will be:
##esm.permdist calpha.permdist MB.general.permdist MB.recessive.permdist MB.dominant.permdist LL.collapse.permdist
write.table( as.data.frame(perms), F,append = TRUE, quote = FALSE, row.names=FALSE, col.names = FALSE )
close(F)

##Release the lock
Unlock(LOCK)
