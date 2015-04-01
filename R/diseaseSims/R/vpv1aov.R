.fillvpv1matrix = function(mm,twoN)
    {
        rv = data.frame((1:twoN)/twoN,NA,NA)
        colnames(rv) <- c("V1","V2","V3")
        for( i in 1:nrow(mm) )
            {
                P = mm[i,1]*twoN
                rv[P,2]=mm[i,2]
                rv[P,3]=mm[i,3]
            }
        ##Courtesy of http://stackoverflow.com/questions/23340150/using-dplyr-window-functions-to-make-trailing-values
        rv %>%
            mutate(dummy = cumsum(0 + !is.na(V2))) %>%
                group_by(dummy,add=TRUE) %>%
                    mutate(filled = nth(V2,1)) %>%
                        mutate(filled2 = nth(V3,1)) %>%
                        ungroup() %>%
                            select(-V2,-V3,-dummy) -> m2
        colnames(m2) <- c("p","rsq","adj.rsq")
        return(m2)
    }

#' Fit linear models to genotypes
#' @param data A list returned by getRiskVariantMatrix
#' @param useSparseM Use SparseM::slm instead of base::lm for the regression.  Requires the SparseM pacakge.
#' @return A data frame with 3 columns: allele frequency, and then estimates of the total variance in the trait explained by markers are <= that frequency.  The estimates make up the last two columns, and are based on the R^2 and adjusted R^2 of the linear model, respectively.
vpv1aov = function(data, useSparseM = FALSE)
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
        data.aov.s = ifelse(USPARSE==FALSE,summary(aov(lm(data$trait ~ ., data = data$genos))),
            summary(aov(SparseM::slm(data$trait ~ ., data=data$genos),data=data$genos)))
        twoN = 2*nrow(data$genos)
        ##Get the counts of risk allele frequencies at each marker
        ##GIVEN THAT THE MARKER WERE USED IN THE REGRESSION
        ROWS=rownames(data.aov.s[[1]])
        alleleCounts = as.integer(colSums(data$genos[,sapply(as.array(ROWS[1:(length(ROWS)-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)]))

        ##Initalize the matrix to return.
        ##The columns will be: mutant allele frequency, r^2, adj. r^2
        p = array( dim = length(unique(alleleCounts)) )
        rsq = array( dim = length(unique(alleleCounts)) )
        adj.rsq = array( dim = length(unique(alleleCounts)) )
        
        sum.sq = data.aov.s[[1]]$'Sum Sq'
        ttl.sum.sq = sum(sum.sq)
        DF = data.aov.s[[1]]$'Df'
        n = length(sum.sq)
        ##Populate the matrix, starting with the rares
        IDX=1
        for( ac in unique(sort(alleleCounts)) )
            {
                ac.sum.sq = sum.sq[which(alleleCounts == ac)]
                p[IDX] = ac/twoN
                ## r^2 due just to mutation at this freq. bin
                rsq[IDX] = sum(ac.sum.sq)/ttl.sum.sq
                if(is.na(rsq[IDX]))
                    {
                        stop("NA found...")
                    }
                ## adj. r^2 due just to mutations at this freq. bin
                adj.rsq[IDX] =  1 - ( (sum.sq[n] + sum(sum.sq[which(alleleCounts != ac)]))/ttl.sum.sq )*(sum(DF)/(DF[n] + length(which(alleleCounts != ac))))
                IDX=IDX+1
            }
        rsq = cumsum(rsq)
        adj.rsq = cumsum(adj.rsq)
        return(.fillvpv1matrix(data.frame(p,rsq,adj.rsq),twoN))
    }
