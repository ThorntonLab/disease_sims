#Takes about 8GB RAM for a c/c panel of 10k individuals
library(SKAT)

##functions
getCCblock=function(con)
    {
        d=readBin(con,"integer",4)
        pos=readBin(con,"numeric",d[3]+d[4])
        genos = matrix(data=0,ncol=d[3]+d[4],nrow=d[1]+d[2])
        for( i in 1:(d[1]+d[2]) )
            {
                nones = readBin(con,"integer",1)
                ones = readBin(con,"integer",nones) + 1
                ntwos = readBin(con,"integer",1)
                twos = readBin(con,"integer",ntwos) + 1
                genos[i,ones] = 1
                genos[i,twos] = 2
            }
        burdens=matrix(readBin(con,"integer",2*(d[1]+d[2])),ncol=2,byrow=TRUE)
        phenos=matrix(readBin(con,"numeric",2*(d[1]+d[2])),ncol=2,byrow=TRUE)
        return(list(pos=pos,genos=genos,burdens=burdens,phenos=phenos))
    }

getSpecificCCblock=function(con,index,recordno)
    {
        z=which(index$V1 == recordno)
        seek(con,index$V2[z])
        return( getCCblock(con) )
    }

mfreq = function(x,ncontrols)
    {
        q = sum(x[1:ncontrols])/(2*ncontrols)
        q = min(q,1-q)
        return (q)
    }


a = commandArgs(trailing = TRUE)

ccindexfile = a[1]
ccfile = a[2]
recno = a[3]
outindex = a[4]
outfile = a[5]
ncontrols = as.integer(a[6])
ncases = as.integer(a[7])
model = a[8]
blocker = a[9]

if( model != "all" & model != "rare" & model != "gwas" )
    {
        warn("Error: invalid model specified")
        q("no")
    }
ccfile.handle = file(ccfile,"rb")

ccindex = read.table(ccindexfile)
ccblock = getSpecificCCblock(ccfile.handle,ccindex,recno)


#generate case/control status indicator variable
status = c(rep(0,ncontrols),rep(1,ncases))

obj = SKAT_Null_Model( status ~ 1,out_type = "D") #dichotomous trait -- see SKAT help(SKAT_Null_Model)

OPIPE=pipe(paste(blocker,recno,outindex,outfile,sep=" "))
if ( model != "all" ) {
                                        #Get minor allele freq (MAF), estimates from controls
    ccblock.maf = apply(ccblock$genos,2,mfreq,ncontrols)
    
    if( model == "rare" ) {
        SKAT.rare = as.numeric( SKAT(ccblock$genos[,which(ccblock.maf < 0.05)], obj)$p.value )
        write(cbind(recno,model,SKAT.rare),OPIPE,ncol=3)
    } else  {
        SKAT.gwas = as.numeric( SKAT(ccblock$genos[,which(ccblock.maf >= 0.05)], obj)$p.value )
        write(cbind(recno,model,SKAT.gwas),OPIPE,ncol=3)
    }
} else {
    SKAT.all = as.numeric( SKAT(ccblock$genos, obj)$p.value )
    write(cbind(recno,model,SKAT.all),OPIPE,ncol=3)
}




