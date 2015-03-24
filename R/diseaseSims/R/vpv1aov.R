.fillvpv1matrix = function(mm,twoN)
    {
        require(dplyr)
        rv = matrix(data=NA,ncol=3,nrow=twoN)
        rv[,1] = (1:twoN)/twoN
        for( i in 1:nrow(mm) )
            {
                P = mm[i,1]*twoN
                rv[P,2]=mm[i,2]
                rv[P,3]=mm[i,3]
            }
        ##Courtesy of http://stackoverflow.com/questions/23340150/using-dplyr-window-functions-to-make-trailing-values
        as.data.frame(rv) %>%
            mutate(dummy = cumsum(0 + !is.na(V2))) %>%
                group_by(dummy,add=TRUE) %>%
                    mutate(filled = nth(V2,1)) %>%
                        mutate(filled2 = nth(V3,1)) %>%
                        ungroup() %>%
                            select(-V2,-V3,-dummy) -> m2
        return(as.matrix(m2))
    }

#' Fit linear models to genotypes
#' @param data A list returned by getRiskVariantMatrix
#' @param ofilename The output file to write to (optional).  File will be gzip-compressed using zlib.
#' @param replicateID An unsigned integer identifying "data".  (optional, only required if ofilename is used)
#' @param append Whether or not to open ofilename in append mode, or to overwrite ofilename
#' @param useSparseM Use SparseM::slm instead of base::lm for the regression.  Requires the SparseM pacakge.
#' @return A matrix with 3 columns: allele frequency, and then estimates of the total variance in the trait explained by markers are <= that frequency.  The estimates make up the last two columns, and are based on the R^2 and adjusted R^2 of the linear model, respectively.
#' @details
#' By default, just the matrix is returned.  If ofilename is specified, the matrix plus a column containing replicateID are also written.  If ofilenname is opened with append = TRUE, then multiple processes are able to write to ofilename, as access is handled using POSIX file locking methods.
vpv1aov = function(data, ofilename , replicateID, append = FALSE, useSparseM = FALSE)
    {
        ##Fit the model and summarize
        ##We coerce the matrix to a data frame so that we can get sum of squares, etc.,
        ##per marker
        USPARSE=useSparseM
        if( USPARSE == TRUE )
            {
                if ( requireNamespace("SparseM",quietly=TRUE) == FALSE )
                    {
                        warning("SparseM not available, using lm() instead")
                        USPARSE=FALSE
                    }
            }
        #data.aov.s = SparseM::slm(data$trait ~ ., data=as.data.frame(data$genos))
        data.aov.s = ifelse(USPARSE==FALSE,summary(aov(lm(data$trait ~ ., data = as.data.frame(data$genos)))),
            summary(aov(SparseM::slm(data$trait ~ ., data=as.data.frame(data$genos)),data=as.data.frame(data$genos))))
        twoN = 2*nrow(data$genos)
        ##Get the counts of risk allele frequencies at each marker
        alleleCounts = as.integer(colSums(data$genos))

        ##Initalize the matrix to return.
        ##The columns will be: mutant allele frequency, r^2, adj. r^2
        mm = matrix(data = NA, ncol = 3, nrow = length(unique(alleleCounts)))

        sum.sq = data.aov.s[[1]]$'Sum Sq'
        ttl.sum.sq = sum(sum.sq)
        DF = data.aov.s$'Df'
        n = length(sum.sq)

        ##Populate the matrix, starting with the rares
        IDX=1
        for( ac in unique(sort(alleleCounts)) )
            {
                ac.sum.sq = sum.sq[which(alleleCounts == ac)]
                mm[IDX,1] = ac/twoN
                ## r^2 due just to mutation at this freq. bin
                mm[IDX,2] = sum(ac.sum.sq)/ttl.sum.sq
                if(is.na(mm[IDX,2]))
                    {
                        stop("NA found...")
                    }
                IDX=IDX+1
            }
        mm[,2] = cumsum(mm[,2])
        if (! missing(ofilename) )
            {
                #.writeVpV1Data(mm,ofilename,replicateID,append);
                .writeVpV1Data(.fillvpv1matrix(mm,twoN),ofilename,replicateID,append)
            }
        return(mm)
    }
