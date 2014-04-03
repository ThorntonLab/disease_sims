readMutsFromPop=function(con,N=1)
  {
    nmuts=readBin(con,"integer",1)
    m=matrix(data=NA,ncol=4,nrow=nmuts,dimnames=list(NULL,c("freq","origin","pos","esize")))
    
    for(i in 1:nmuts)
      {
        ID=readBin(f,"integer",1)
        freq=readBin(f,"integer",1)/N
        gen=readBin(f,"integer",1)
        neutral=readBin(f,"integer",1,size=1)
        pos=readBin(f,"numeric",1)
        esize=readBin(f,"numeric",1)
        label=readBin(f,integer(),1,size=1)
        
        m[i,1]=freq
        m[i,2]=gen
        m[i,3]=pos
        m[i,4]=esize
      }
    return (as.data.frame(m))
  }
