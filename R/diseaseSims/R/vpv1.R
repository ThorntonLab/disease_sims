#' Calculate additive variation explained as a function of mutation frequency
#' @param data A list that is the return value from getRiskVariantMatrix
#' @param ofilename Output file name.  (Optional)
#' @param replicateID A unique identifier for this data set (unsigned integers, please!!)
#' @param append Append to output file?
#' @param maxfreq Maximum frequency to work up to
#' @return fitted.vals, which is a list of fitted values for each marker
#' @return mm, which is a 3-column matrix: p, variance explained by all markers with frequency <= p, based on R^2, and then variance explained by all markers with frequency <= p, based on adjusted R^2
#' @details
#' Output to the access file is managed via a POSIX file lock.  This lock scheme allows you to analyze multiple data sets
#' using your cluster, and then write all the output to a single file.  To make this work, though, be sure that append == TRUE,
#' otherwise you'll lock the file and then over-write it with each process!
vpv1 = function( data, ofilename , replicateID, append = FALSE, maxfreq = 1.0 )
    {
        if ( maxfreq < 0 | maxfreq > 1 )
            stop("Error, maxfreq must be >= 0 and <= 1")
        twoN = 2*nrow(data$genos)
        data.genos.S = as.integer(colSums(data$genos))
        MVAL = min(as.integer(maxfreq*twoN),twoN)
        mm = matrix(data=NA,ncol=3,nrow=MVAL)
        LLold=array()
        for( P in 1:MVAL )
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
                        lm_i.s = summary(RcppEigen::fastLm( data$trait ~ data$genos[,LL]))
                        #lm_i.s = summary(lm( data$trait ~ data$genos[,LL]))
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
        ll = lm( data$trait ~ data$genos )
        ll.fittedVals = fitted(ll)
        if (! missing(ofilename) )
            {
                .writeVpV1Data(mm,ofilename,replicateID,append);
            }
        return( list(fitted.vals=ll.fittedVals,vexp=mm) )
    }

#' Some stuff
vp = function( data, ofilename , replicateID, freqrange, append = FALSE )
    {
        twoN = 2*nrow(data$genos)
        data.genos.S = as.integer(colSums(data$genos))
        mm = matrix(data=NA,ncol=3,nrow=length(freqrange))
        LLold=array()
        IDX=1
        for( p in freqrange )
            {
                P = as.integer( p*twoN )
                LL=which(data.genos.S <= P)
                mm[IDX,1] = P/twoN
                LENDIFF = ifelse(length(LLold)!=length(LL),TRUE,FALSE)
                DATADIFF = ifelse( LENDIFF == TRUE, TRUE, FALSE )
                if( LENDIFF == FALSE )
                    {
                        DATADIFF = ifelse( length(which(LLold == LL)) != length(LL), TRUE,FALSE )
                    }
                if( LENDIFF | DATADIFF )
                    {
                        lm_i.s = summary(RcppEigen::fastLm( data$trait ~ data$genos[,LL]))
                        #lm_i.s = summary(lm( data$trait ~ data$genos[,LL]))
                        mm[IDX,2] = lm_i.s$r.squared
                        mm[IDX,3] = lm_i.s$adj.r.squared
                    }
                else
                    {
                        mm[IDX,2] = mm[IDX-1,2]
                        mm[IDX,3] = mm[IDX-1,3]
                    }
                LLold=LL
                IDX=IDX+1
            }
        if (! missing(ofilename) )
            {
                .writeVpV1Data(mm,ofilename,replicateID,append);
            }
        return(mm)
    }
