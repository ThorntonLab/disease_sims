library(diseaseSims)
library(buRden)
library(data.table)

n = commandArgs(trailing = TRUE)

## Replicate ID
repID = as.integer(n[1])
CCIDX=n[2]
CCFILE=n[3]
PERMFILE=n[4]
## output file name
OFILENAME = n[5]

perms = fread(paste("gunzip -c ",PERMFILE))

setnames(perms,c("esm.permdist","calpha.permdist","MB.general.permdist", "MB.recessive.permdist", "MB.dominant.permdist", "LL.collapse.permdist"))

##Get the observed data
ccidx = read.table(CCIDX,colClasses=c("integer","numeric"))
ccdata = getCCblock(CCFILE,procIndex(ccidx,2,repID))
status = c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))

##Minor allele frequencies in controls
MAF = colSums(ccdata$genos[1:ccdata$ncontrols,])/(2*ccdata$ncontrols)
##number of markers for the ESM test
K=50

LLCMAF = 0.05

##Get the observe values
obs = allBurdenStats(ccdata$genos[,which(MAF < 0.05)],status,K,LLCMAF)
NPERMS = nrow(perms)
##Now, get one tailed Pr(X >= x), Z-score p-values, and store in data.frame for easy output
esm.z = (obs$esm.stat - mean(perms$esm.permdist) )/sd(perms$esm.permdist)
calpha.z = (obs$calpha.stat - mean(perms$calpha.permdist))/sd(perms$calpha.permdist)
MB.g.z = (obs$MB.general.stat-mean(perms$MB.general.permdist))/sd(perms$MB.general.permdist)
MB.r.z = (obs$MB.recessive.stat-mean(perms$MB.recessive.permdist))/sd(perms$MB.recessive.permdist)
MB.d.z = (obs$MB.dominant.stat-mean(perms$MB.dominant.permdist))/sd(perms$MB.dominant.permdist)
LL.z = (obs$LL.collapse.stat-mean(perms$LL.collapse.permdist))/sd(perms$LL.collapse.permdist)
PVALS = data.frame(esm.perm.p = length(which( perms$esm.permdist >= obs$esm.stat ))/NPERMS,
    esm.z.p = 2*pnorm(abs(esm.z),lower=F),
    calpha.perm.p = length(which(perms$calpha.permdist >= obs$calpha.stat))/NPERMS,
    calpha.z.p = 2*pnorm(abs(calpha.z),lower=F),
    MB.general.p = length(which(perms$MB.general.permdist >= obs$MB.general.stat))/NPERMS,
    MB.general.z.p = 2*pnorm(abs(MB.g.z),lower=F),
    MB.recessive.p = length(which(perms$MB.recessive.permdist >= obs$MB.recessive.stat))/NPERMS,
    MB.recessive.z.p = 2*pnorm(abs(MB.r.z),lower=F),
    MB.dominant.p = length(which(perms$MB.dominant.permdist >= obs$MB.dominant.stat))/NPERMS,
    MB.dominant.z.p = 2*pnorm(abs(MB.d.z),lower=F),
    LL.collapse.p = length(which(perms$LL.collapse.permdist >= obs$LL.collapse.stat))/NPERMS,
    LL.z.p = 2*pnorm(abs(LL.z),lower=F))

LOCK = Lock(OFILENAME)

OF = gzfile(OFILENAME,"w")

write.table(PVALS,OF,col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

close(OF)

Unlock(LOCK)


