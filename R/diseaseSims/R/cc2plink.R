#' reformat a case/control block to ped/map blocks
cc2plink=function(ccblock,causalOnly = FALSE) {
    ped = data.frame( cbind(familyID = rep(1,nrow(ccblock$genos)),
        indID = 1:nrow(ccblock$genos),
        paternalID = rep(0,nrow(ccblock$genos)),
        maternalID = rep(0,nrow(ccblock$genos)),
        sex = rep(0,nrow(ccblock$genos)),
        phenotype = c(rep(1,ccblock$ncontrols),rep(2,ccblock$ncases)),
        .reformatCCgenos(ccblock$genos)) )
    map = data.frame("chrom"=rep(1,length(ccblock$pos)),
        "snpID"=paste("rs",1:length(ccblock$pos),sep=""),
        "genetic"=ccblock$pos,
        "physical"=ccblock$pos)
    return(list(ped=ped,map=map))
}
