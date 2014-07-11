#note: this takes ~1GB RAM for 5000 cases & 5000 controls

#gives case/control p-values for single-marker test under an additive model
n=commandArgs(trailing=T)

indexfile=n[1]
anovaindex=n[2]
effectfile=n[3]
recordno=as.integer(n[4])  #"type cast" needed here, else passing to writePVblock returns undecipherable error :)
anovafile=n[5]
ncontrols=as.integer(n[6])
ncases=as.integer(n[7])
N=as.integer(n[8])
outfile=n[9]
outindex=n[10]
blocker=n[11]

index=read.table(indexfile)
aindex=read.table(anovaindex)

esizes = getEsizes(effectfile,procIndex(index,2,recordno))
ccdata = getCCblock(anovafile,procIndex(aindex,2,recordno))

FAILFIND=array()
FAILNO=0
for( i in 1:length(ccdata$pos) )
    {
        if( length(which(esizes[,1] == ccdata$pos[i])) == 0 )
            {
                FAILNO=FAILNO+1
                FAILFIND[FAILNO]=ccdata$pos[i]
            }
    }

if(FAILNO>0)
{
    stop(paste("Fatal error: ", FAILNO," mutations in ",anovafile,
               " were not found in ",effectfile))
}

status=c(rep(0,ccdata$ncontrols),rep(1,ccdata$ncases))
pvblock=makePVblock(ccdata,esizes,status)
pvblock[,"popfreq"]=pvblock[,"popfreq"]/(2*N)

#OLD CODE:
#options(warn=-1)
#f=pipe(paste(blocker,recordno,outindex,outfile,as.integer(nrow(pvblock)),sep=" "))
#write.table(pvblock,f,row.names=F,col.names=T,quote=F,append=T)

#Uses Rcpp to manage POSIX file locking from within R:
writePVblock(outfile,outindex,recordno,pvblock)


