n=commandArgs(trailing=TRUE)

indexfile=n[1]
pvalfile=n[2]

x=read.table(indexfile)

f=file(pvalfile,"r")

rejrate=0
for(i in 1:nrow(x))
    {
        seek(f,x$V2[i])
        data=read.table(f,header=T,nrows=x$V3[i])
        if ( length(which(data$score >= 8)) > 0 )
            {
                rejrate=rejrate+1
            }
    }
print(rejrate/nrow(x))
