#' LD summaries of simulated case/control panels
#' @param ccblock A case/control panel retured by diseaseSims::getCCblock
#' @param pvblock Single-marker association results from diseaseSims::makePVblock
#' @param sig.threshold The -log10(p-value) beyond which a single-marker test is considered "significant"
#' @param minrsq Only include causative markers in the return value if their LD with the significant marker is >= minrsq
#' @return A data frame summarizing LD beteween significant markers and causative markers.  Included are the population frequencies of the mutations and the effect sizes of the causative mutations.  The data frame contains: sp = position of the significant marker.  cp = position of causative marker.  sfreq = frequency of significant mutation in general population.  cfreq = frequency of causative mutation in general population. ce = effect size of causative mutation.  rsq = squared correlation coefficient between significant and causative marker genotype.
ccLD=function(ccblock,pvblock,sig.threshold = 8,minrsq = 0 )
    {
        if( ncol(ccblock$genos) != nrow(pvblock) )
           {
               stop("ccLD error: nrow(ccblock$genos) != nrow(pvblock)")
           }
        if( minrsq < 0 || minrsq > 1 )
            {
                stop("ccLD error: minrsq must be <= 0 and <= 1")
            }
        sigpos = which(pvblock$score >= sig.threshold)
        causative = which( pvblock$esize != 0 )
        if( length(causative) != ccblock$causative )
            {
                stop("ccLD error: number of causative mutation in pvblock != number in ccbloc$genos");
            }
        #Arrays to store return value
        sp=array()
        cps=array()
        sfreq=array()
        cfreq=array()
        ce=array()
        LD=array()
        DUMMY=1
        if ( length(sigpos) > 0 && length(causative) > 0 )
            {
                for( sigp in 1:length(sigpos) )
                    {
                        for(cp in 1:length(causative))
                            {
                                #Genotype correlation coefficient = the standard r^2
                                rsq = cor( ccblock$genos[,sigpos[sigp]],ccblock$genos[,causative[cp]] )^2
                                if( rsq >= minrsq )
                                    {
                                        sp[DUMMY]=pvblock$pos[sigpos[sigp]]
                                        cps[DUMMY]=pvblock$pos[causative[cp]]
                                        sfreq[DUMMY]=pvblock$popfreq[sigpos[sigp]]
                                        cfreq[DUMMY]=pvblock$popfreq[causative[cp]]
                                        ce[DUMMY]=pvblock$esize[causative[cp]]
                                         LD[DUMMY]=rsq
                                        DUMMY=DUMMY+1
                                    }
                            }
                    }
            }
        return( data.frame(sigpos = sp,
                           cpos = cps,
                           sigfreq = sfreq,
                           cfreq = cfreq,
                           esize = ce,
                           rsq = LD ) )
    }
