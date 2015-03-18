#' Calculate additive variation explained as a function of mutation frequency
#' @param data A list that is the return value from getRiskVariantMatrix
#' @return fitted.vals, which is a list of fitted values for each marker
#' @return mm, which is a 3-column matrix: p, variance explained by all markers with frequency <= p, based on R^2, and then variance explained by all markers with frequency <= p, based on adjusted R^2
#' @details
#' This is a general method that goes beyond the simple formula presented in the literature, which is based on additve models
vpv1 = function( data )
    {
        twoN = 2*nrow(data$genos)
        data.genos.S = as.integer(colSums(data$genos))
        mm = matrix(data=NA,ncol=3,nrow=twoN)
        LLold=array()
        for( P in 1:twoN )
            {
                LL=which(data.genos.S <= P)
                mm[P,1] = P/twoN
                LENDIFF = ifelse(length(LLold)!=length(LL),TRUE,FALSE)
                DATADIFF = ifelse( LENDIFF == TRUE, TRUE, FALSE )
                if( LENDIFF == FALSE )
                    {
                        DATADIFF = ifelse( length(which(LLold == LL)) != length(LL), TRUE,FALSE )
                    }
                if( LENDIFF | DATADIFF )
                    {
                        lm_i = lm( data$G ~ data$genos[,LL])
                        lm_i.s = summary(lm_i)
                        mm[P,2] = lm_i.s$r.squared
                        mm[P,3] = lm_i.s$adj.r.squared
                    }
                else
                    {
                        mm[P,2] = mm[P-1,2]
                        mm[P,3] = mm[P-1,3]
                    }
                LLold=LL
            }
        ll = lm( data$G ~ data$genos )
        ll.fittedVals = fitted(ll)
        return( list(fitted.vals=ll.fittedVals,vexp=mm) )
    }
