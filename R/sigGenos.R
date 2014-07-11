#Get genotypes at significantly-associated markers
library(diseaseSims)

n=commandArgs(trailing=TRUE)

smpvfile=n[1]
smindexfile=n[2]
ccfile=n[3]
ccindexfile=n[4]
ccpeoplefile=n[4]

smindex=read.table(smindexfile)
ccindex=read.table(ccindexfile)

smfile=file(smpvfile,"r")
for( i in 1:nrow(smindex) )
{
    recordno=as.integer(smindex$V1[i]);
    data=read.table(smfile,header=T,nrow=smindex$V3[i])
    ccblock=getCCblock(ccfile,procIndex(ccindex,2,recordno))
}
