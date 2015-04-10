##Calculate true heritability, etc., in the case/control panel

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

CONTRANGE = 1:ccdata$ncontrols
CASERANGE = (ccdata$ncontrols+1):nrow(ccdata$genos)

print(range(CONTRANGE))
print(range(CASERANGE))
cch2 = list( total.vg = var(ccdata$phenos[,1]),
    total.ve = var(ccdata$phenos[,2]),
    controls.vg = var(ccdata$phenos[CONTRANGE,1]),
    controls.ve = var(ccdata$phenos[CONTRANGE,2]),
    cases.vg = var(ccdata$phenos[CASERANGE,1]),
    cases.ve = var(ccdata$phenos[CASERANGE,2]),
    vg.total = var(ccdata$phenos[,1]),
    ve.total = var(ccdata$phenos[,2]),
    v.pheno = var(ccdata$phenos[,1]+ccdata$phenos[,2])
    )

cch2["total.h2"] = cch2$total.vg/(cch2$total.vg +cch2$total.ve)
cch2["controls.h2"] = cch2$controls.vg/(cch2$controls.vg+cch2$controls.ve)
cch2["cases.h2"] = cch2$cases.vg/(cch2$cases.vg+cch2$cases.ve)

print(cch2)
print(cch2$vg.total/cch2$v.pheno)



