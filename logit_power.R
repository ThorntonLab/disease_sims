n=commandArgs(trailing=TRUE)

indexfile=n[1]
pvalfile=n[2]

x=read.table(indexfile)

f=file(pvalfile,"r")

rejrate=0
rejrateGWAS=0
for(i in 1:nrow(x))
    {
        seek(f,x$V2[i])
        data=read.table(f,header=T,nrows=x$V3[i])
        if ( length(which(data$score >= 8)) > 0 )
            {
                rejrate=rejrate+1
            }
        z=which(data$popfreq >= 0.05 & data$popfreq <= 0.95)
        if ( length(which(data$score[z] >= 8)) > 0 )
            {
                rejrateGWAS=rejrateGWAS+1
            }
    }
print(paste(rejrate/nrow(x),rejrateGWAS/nrow(x)))
