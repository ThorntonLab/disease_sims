#' Returns the indexes of the top n hits
#' @param x A vector
#' @param n The number of top hits that you wish to return
#' @return which(rank(1/abs(x)) <= min(n,length(x))) )
bestHits = function(x,n)
    {
        return( which(rank(1/abs(x)) <= min(n,length(x))) )
    }

#' LD summaries of simulated case/control panels
#' @param ccblock A case/control panel retured by diseaseSims::getCCblock
#' @param pvblock Single-marker association results from diseaseSims::makePVblock
#' @param sig.threshold The -log10(p-value) beyond which a single-marker test is considered "significant"
#' @param minrsq Only include causative markers in the return value if their LD with the significant marker is >= minrsq
#' @param maxvals If > 0, return the top maxvals records per marker.  In other words, return the maxvals highest associations per significant marker.  maxvals = 1 returns the top hit, 2 returns the top 2, etc.
#' @return A data frame summarizing LD beteween significant markers and causative markers.  Included are the population frequencies of the mutations and the effect sizes of the causative mutations.  The data frame contains: sp = position of the significant marker.  cp = position of causative marker.  sfreq = frequency of significant mutation in general population.  cfreq = frequency of causative mutation in general population. sesize = effect size of significant mutation. cesize = effect size of causative mutation.  rsq = squared correlation coefficient between significant and causative marker genotype.
ccLD=function(ccblock,pvblock,sig.threshold = 8,minrsq=0,maxvals=-1)
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
        se=array()
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
                                        se[DUMMY]=pvblock$esize[sigpos[sigp]]
                                        ce[DUMMY]=pvblock$esize[causative[cp]]
                                        LD[DUMMY]=rsq
                                        DUMMY=DUMMY+1
                                    }
                            }
                    }
            }
        
        if( maxvals == -1 )
            {
                rv = data.frame(sigpos = sp,
                    cpos = cps,
                    sigfreq = sfreq,
                    cfreq = cfreq,
                    sesize = se,
                    cesize = ce,
                    rsq = LD )
                return (rv);
            }
        else
            {
                rv = data.frame(sigpos = sp,
                    cpos = cps,
                    sigfreq = sfreq,
                    cfreq = cfreq,
                    sesize = se,
                    cesize = ce,
                    rsq = LD,
                    mx = rep(maxvals,length(sp)))
                
                rv.mx = plyr::ddply( rv,.(sigpos),
                    summarise,
                    cpos = cpos[ bestHits(rsq,mx) ],
                    sigfreqs = sigfreq[ bestHits(rsq,mx) ],
                    cfreq = cfreq[ bestHits(rsq,mx) ],
                    sesize = sesize[ bestHits(rsq,mx) ],
                    cesize = cesize[ bestHits(rsq,mx) ],
                    rsq = rsq[ bestHits(rsq,mx) ] )
                return(rv.mx)
            }
    }
