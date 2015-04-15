#' reformat a case/control block to ped/map blocks
cc2plink=function(ccblock,causalOnly = FALSE,phenos = FALSE) {
    P = array();
    S=array()
    if( phenos == FALSE ) {
        P=c(rep(1,ccblock$ncontrols),rep(2,ccblock$ncases))
        S = rep(0,nrow(ccblock$genos))
    } else {
        P = rowSums(ccblock$phenos)
        S = rep(1,nrow(ccblock$genos))
    }
    ped=data.frame()
    map=data.frame()
    if(causalOnly==FALSE)
        {
    ped = data.frame( cbind(familyID = rep(1,nrow(ccblock$genos)),
        indID = 1:nrow(ccblock$genos),
        paternalID = rep(0,nrow(ccblock$genos)),
        maternalID = rep(0,nrow(ccblock$genos)),
        sex = S,
        phenotype = P,
        .reformatCCgenos(ccblock$genos)) )
    map = data.frame("chrom"=rep(1,length(ccblock$pos)),
        "snpID"=paste("rs",1:length(ccblock$pos),sep=""),
        "genetic"=ccblock$pos,
        "physical"=ccblock$pos)
} else
    {
    Z=(ccblock$neutral+1):ncol(ccblock$genos)
    ped = data.frame( cbind(familyID = rep(1,nrow(ccblock$genos)),
        indID = 1:nrow(ccblock$genos),
        paternalID = rep(0,nrow(ccblock$genos)),
        maternalID = rep(0,nrow(ccblock$genos)),
        sex = S,
        phenotype = P,
        .reformatCCgenos(ccblock$genos[,Z])))
    map = data.frame("chrom"=rep(1,length(ccblock$pos[Z])),
        "snpID"=paste("rs",1:length(ccblock$pos[Z]),sep=""),
        "genetic"=ccblock$pos[Z],
        "physical"=ccblock$pos[Z])
}
    return(list(ped=ped,map=map))
}
