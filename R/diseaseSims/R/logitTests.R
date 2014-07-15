#' Association testing under an additive model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.additive=function(x,status)
{
    if(length(x) != length(status))
        {
            stop("logit.additive: length(x) != length(status)")
            return;
        }
    GLM=glm( status ~ x, family=binomial("logit") )
    return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
}

#' Association testing under a recessive model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @note The test works by recoding 1 genotypes as 0 (e.g., one copy of minor allele is the same as none)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.recessive=function(x,status)
    {
        if(length(x) != length(status))
            {
                stop("logit.additive: length(x) != length(status)")
                return;
            }
        x[which(x==1)]=0
        GLM=glm( status ~ x, family=binomial("logit") )
        return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
    }

#' Association testing under a dominant model
#' @param x An array of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal length(x)
#' @note The test works by recoding 2 genotypes as 1 (e.g., one copy of minor allele is the same as two)
#' @return The p-value of a logistic regression of case/control status onto genotype.
logit.dominant=function(x,status)
    {
        if(length(x) != length(status))
            {
                stop("logit.dominant: length(x) != length(status)")
                return;
            }
        x[which(x==2)]=1
        GLM=glm( status ~ x, family=binomial("logit") )
        return( as.numeric( summary(GLM)$coefficients[,4][2] ) )
    }

#' Obtain single-marker p-values for case/control data
#' @param genos A matrix of genotypes, coded as number of copies of minor allele
#' @param status An array of 0 = control, 1 = case.  length(status) must equal nrow(x)
#' @param model One of "additive","recessive", or "dominant"
#' @return An array of p-values of logistic regressions of case/control status onto genotype.  The order of the p-values corresponds to the column order in x.
ccpvals = function(genos,status,model="additive")
    {
        if( length(status) != nrow(genos) )
            {
                stop("ccpvals.additive: nrow(x) != length(status)")
                return;
            }
        if( model == "additive" )
            {
                return(apply(genos,2,logit.additive,status))
            }
        else if ( model == "recessive" )
            {
                return(apply(genos,2,logit.recessive,status))
            }
        else if (model == "dominant")
            {
                return(apply(genos,2,logit.dominant,status))
            }
        else
            {
                stop(paste("ccpvals: model",model,"is not valid. Model must be one of additive, dominant, or recessive"))
            }
        return();
    }

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
#' #Make a random matrix of genotypes for 500 controls, 500 cases, 100 markers
#' set.seed(100)
#' ncases=500
#' ncontrols=500
#' nmarkers=100
#' genos=matrix(data=replicate(nmarkers,rbinom(ncases+ncontrols,2,0.5)),nrow=ncases+ncontrols,ncol=nmarkers,byrow=FALSE)
#' #Make a fake esizes block:
#' pos=runif(nmarkers,0,1)
#' esize=c(rep(0,80),rexp(20))
#' count=as.integer(runif(nmarkers,1,1e3))
#' ages=count
#' esizes=list(pos=pos,esize=esize,count=count,age=ages)
#' #This is an incomplete "ccblock", but good enough to run through this function:
#' ccblock=list(genos=genos,pos=pos,ncases=ncases,ncontrols=ncontrols)
#' ccblock.pvblock = makePVblock(ccblock,esizes,c(rep(0,ncontrols),rep(1,ncases)))
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
        output[,"score"]=-log10(ccpvals(ccdata$genos,status,model))
        return(as.data.frame(output))
    }
