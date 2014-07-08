#note: this takes ~1GB RAM for 5000 cases & 5000 controls

#gives case/control p-values for single-marker test under an additive model

#getPheno=function(con,N)
#    {
#      nphenos=readBin(con,"integer",1)
#      pheno=matrix( readBin(con,"numeric",2*N), ncol=2, byrow=TRUE,dimnames=list(NULL,c("G","E")))
#      haps=matrix( readBin(con,"integer",2*N), ncol=2, byrow=TRUE,dimnames=list(NULL,c("h1","h2")))
#      return( list(pheno=pheno,haps=haps) )
#    }

#getSpecificPheno=function(con,index,recordno,N)
#    {
#        z=which(index$V1 == recordno)
#        seek(con,index$V3[z])
#        return (getPheno(con,N))
#    }

#getEsizes=function(con)
#    {
#        nm=readBin(con,"integer",1)
#        return( as.data.frame(matrix( readBin(con,"numeric",4*nm), ncol=4,byrow=TRUE,
#                                     dimnames=list(NULL,c("pos","esize","count","age")))) )
#    }

getEoffset = function(indexfile,recordno)
    {
        idx=read.table(indexfile);
        z=which(idx$V1==recordno)
        #print(z)
        if( z==1 )
            {
                return (0);
            }
        ispos = if (length(which(diff(idx$V2[1:(z-1)])>0)) == z - 2) TRUE else FALSE
        if( ispos )
            {
                return (idx$V2[z])
            }
        else
            {
                return( sum(idx$V2[1:(z-1)]) )
            }
    }

getSpecificEsizes=function(filename,index,recordno)
    {
        offset = getEoffset(index,recordno)
        print(paste("e:",offset))
        return ( getEsizes( filename,offset ) )
    }

#getCCblock=function(con)
#    {
#        d=readBin(con,"integer",4)
#        pos=readBin(con,"numeric",d[3]+d[4])
#        genos = matrix(data=0,ncol=d[3]+d[4],nrow=d[1]+d[2])
#        for( i in 1:(d[1]+d[2]) )
#            {
#                nones = readBin(con,"integer",1)
#                ones = readBin(con,"integer",nones) + 1
#                ntwos = readBin(con,"integer",1)
#                twos = readBin(con,"integer",ntwos) + 1
#                genos[i,ones] = 1
#                genos[i,twos] = 2
#            }
#        burdens=matrix(readBin(con,"integer",2*(d[1]+d[2])),ncol=2,byrow=TRUE)
#        phenos=matrix(readBin(con,"numeric",2*(d[1]+d[2])),ncol=2,byrow=TRUE)
#        return(list(pos=pos,genos=genos,burdens=burdens,phenos=phenos))
#    }

getCCoffset = function(indexfile,recordno)
    {
        #cheat: the relevant column is the second
        return (getEoffset(indexfile,recordno))
    }
    
getSpecificCCblock=function(filename,index,recordno)
    {
        offset = getCCoffset(index,recordno)
        print(paste("cc:",offset))
        return( getCCblock(filename,offset) )
    }

dologit=function(x,ncontrols,ncases)
{
    status = c(array(dim=ncontrols,0),array(dim=ncases,1))
    GLM=glm( status ~ x, family=binomial("logit") )
    return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
}

getPvalsCCstatus = function(genos,ncontrols,ncases)
  #logistic regression of case control status on
  #genotype as an additive model, as per Goldstein's simulation study
    {
        return(apply(genos,2,dologit,ncontrols,ncases))
    }

## readMutsFromPop=function(con,N=1)
##   {
##     nmuts=readBin(con,"integer",1)
##     m=matrix(data=NA,ncol=4,nrow=nmuts,dimnames=list(NULL,c("freq","origin","pos","esize")))
    
##     for(i in 1:nmuts)
##       {
##         ID=readBin(f,"integer",1)
##         freq=readBin(f,"integer",1)/N
##         gen=readBin(f,"integer",1)
##         neutral=readBin(f,"integer",1,size=1)
##         pos=readBin(f,"numeric",1)
##         esize=readBin(f,"numeric",1)
##         label=readBin(f,integer(),1,size=1)
        
##         m[i,1]=freq
##         m[i,2]=gen
##         m[i,3]=pos
##         m[i,4]=esize
##       }
##     return (as.data.frame(m))
##   }

## readSpecificMutsFromPop=function(con,N,index,recordno)
##   {
##     z=which(index$V1 == recordno)
##     seek(con,index$V2[z])
##     return (readMutsFromPop(con,N))
##   }

makePVblock = function( ccdata, esizes, ncontrols,ncases )
    {
        output=matrix(data=NA,ncol=6,nrow=ncol(ccdata$genos),
            dimnames=list(NULL,(c("pos","esize","mfcontrols","mfcases","popfreq","score"))))
        output[,"pos"]=ccdata$pos
        output[,"mfcontrols"]=colSums(ccdata$genos[1:ncontrols,])/(2*ncases)
        output[,"mfcases"]=colSums(ccdata$genos[(ncontrols+1):nrow(ccdata$genos),])/(2*ncases)
        for( r in 1:nrow(output) )
            {
                z=which(as.numeric(esizes$pos) == ccdata$pos[r]);
                #if(length(z)==0)
                #    {
                        #mutation is neutral, get its frequency from the population
                #        output[r,"esize"]=0
                #        output[r,"popfreq"]=NA
                #    }
                #else
                #    {
                output[r,"esize"]=as.numeric(esizes$esize[z])
                output[r,"popfreq"]=as.numeric(esizes$count[z])
                #    }
            }
        output[,"score"]=-log10(getPvalsCCstatus(ccdata$genos, ncontrols,ncases))
        return(output)
    }

n=commandArgs(trailing=T)

indexfile=n[1]
anovaindex=n[2]
effectfile=n[3]
recordno=n[4]
anovafile=n[5]
ncontrols=as.integer(n[6])
ncases=as.integer(n[7])
N=as.integer(n[8])
outfile=n[9]
outindex=n[10]
blocker=n[11]
#index=read.table(indexfile)
#aindex=read.table(anovaindex)

#f=file(effectfile,"rb")
#esizes=readSpecificMutsFromPop(f,2*N,index,recordno)
#esizes=if(grep("gz",effectfile)) getSpecificEsizesGZ(f,index,recordno) else getSpecificEsizes(f,index,recordno)
esizes=getSpecificEsizes(effectfile,indexfile,recordno)
print(nrow(esizes))
#close(f)
#f=file(anovafile,"rb")
ccdata=getSpecificCCblock(anovafile,anovaindex,recordno)

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

pvblock=makePVblock(ccdata,esizes,ncontrols,ncases)
pvblock[,"popfreq"]=pvblock[,"popfreq"]/(2*N)

options(warn=-1)
f=pipe(paste(blocker,recordno,outindex,outfile,as.integer(nrow(pvblock)),sep=" "))
write.table(pvblock,f,row.names=F,col.names=T,quote=F,append=T)


