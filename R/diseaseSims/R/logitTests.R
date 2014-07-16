#' Perform single-marker association test on simulated case/control data
#' @param ccdata A case/control data set returned from \link[diseaseSims]{getCCblock}
#' @param esizes A matrix summarizing mutation effect sizes returned from \link[diseaseSims]{getEsizes}
#' @param status An array specifying discrete phenotypes (0=control,1=case).  Length status must equal nrow(ccdata$geno)
#' @param model Model for the association test.  Must be one of additive, recessive, or domimant
#' @return A data frame containing the following columns:
#' pos = the position of the mutation,
#' esize = the effect size of the mutation,
#' mfcontrols = minor allele frequency in controls,
#' mfcases = minor allele frequency in cases,
#' popfreq = DERIVED allele frequency in entire population (NOT minor allele frequency!!!)
#' score = -log10(logistic regression p-value)
#' @examples
#' data(ccblock.dsdata)
#' data(esizes.dsdata)
#' pvblock = makePVblock(ccblock.dsdata,esizes.dsdata,c(rep(0,ccblock.dsdata$ncontrols),rep(1,ccblock.dsdata$ncases)))
makePVblock = function( ccdata,
    esizes,
    status,
    model = "additive")
    {
        ncases = length(which(status==0))
        ncontrols = length(which(status==1))
        output=matrix(data=NA,ncol=6,nrow=ncol(ccdata$genos),
            dimnames=list(NULL,(c("pos","esize","mfcontrols","mfcases","popfreq","score"))))
        output[,"pos"]=ccdata$pos
        output[,"mfcontrols"]=colSums(ccdata$genos[which(status==0),])/(2*ncontrols)
        output[,"mfcases"]=colSums(ccdata$genos[which(status==1),])/(2*ncases)
        for( r in 1:nrow(output) )
            {
                z=which(as.numeric(esizes$pos) == ccdata$pos[r]);
                output[r,"esize"]=as.numeric(esizes$esize[z])
                output[r,"popfreq"]=as.numeric(esizes$count[z])
            }
        output[,"score"]=-log10(buRden::ccpvals(ccdata$genos,status,model))
        return(as.data.frame(output))
    }
