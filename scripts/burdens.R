library(diseaseSims)
library(buRden)

ccidx = read.table("ccindex.txt",colClasses=c("integer","numeric"))
ccdata = getCCblock("ccfile.bin.gz",procIndex(ccidx,2,0))
status = c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))


##For fun, let's only run tests on markers with MAF < 0.05 in controls
MAF = colSums(ccdata$genos[1:ccdata$ncontrols,])/(2*ccdata$ncontrols)

##Let's permute all available tests in the buRden package 10 times
NPERMS=10
all.p = allBurdenStats.p.perm(ccdata$genos[,which(MAF<0.05)],status,NPERMS,50,0.025)

##Print out the z-scores:
print(paste("ESM cAlpha MBg MBr MBd LLc"))
print(paste(all.p$esm.z.value,
            all.p$calpha.z.value,
            all.p$MB.general.z.value,
            all.p$MB.recessive.z.value,
            all.p$MB.dominant.z.value,
            all.p$LL.collapse.z.value))


