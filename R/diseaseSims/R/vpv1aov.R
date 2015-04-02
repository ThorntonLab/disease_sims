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

.dobigaov = function(data,chunksize = 5000)
    {
        if( requireNamespace("biglm",quietly=TRUE) == FALSE )
            {
                stop("cannot proceed: biglm namespace not found")
            }
        FINAL = as.numeric(nrow(data))
        IIa=1
        IIb=min(chunksize,FINAL)
        ##Make the model expression
        FORMULA <- paste(colnames(data)[1], "~", paste(colnames(data)[-1], collapse=" + "))
        BIGGIE = biglm::biglm( as.formula(FORMULA) ,data = data[IIa:IIb,] )
        IIa = IIb + 1
        IIb = min(IIa + chunksize - 1, FINAL)
        while( as.numeric(IIa) < FINAL )
            {
                update( BIGGIE, moredata=data[IIa:IIb,] )
                IIa = IIb + 1
                IIb = min(IIa + chunksize - 1, FINAL)
            }
        RV =summary(aov(BIGGIE,data=data))
        return(RV)
    }

#' Fit linear models to genotypes
#' @param data A list returned by getRiskVariantMatrix
#' @param method Method to use for model fitting.  Must be one of "biglm","slm", or "lm", corresponding to the biglm, SparseM, and base packages, respectively.
#' @param chunksize If the data set is massive, we'll use biglm::biglm for the regression, breaking data up into sets of chunksize rows
#' @return A data frame with 3 columns: allele frequency, and then estimates of the total variance in the trait explained by markers are <= that frequency.  The estimates make up the last two columns, and are based on the R^2 and adjusted R^2 of the linear model, respectively.
vpv1aov = function(data, method = "lm", chunksize=5000)
    {
        ##Check if data$geno's dimensions are potentially "too big":
        BIG = ifelse( as.numeric(ncol(data$genos))*as.numeric(nrow(data$genos)) >= as.numeric(.Machine$integer.max),
            TRUE,
            ifelse( method == "biglm", TRUE, FALSE) )

        if( BIG == TRUE & method != "biglm" )
            {
                warning("data are very big, will use biglm")
            }
        ##Fit the model and summarize
        ##We coerce the matrix to a data frame so that we can get sum of squares, etc.,
        ##per marker
        USPARSE = ifelse( method == "slm",
            ifelse( BIG == FALSE, TRUE, FALSE ),
            FALSE )
        if( USPARSE == TRUE )
            {
                if ( requireNamespace("SparseM",quietly=TRUE) == FALSE )
                    {
                        warning("SparseM not available, using lm() instead")
                        USPARSE=FALSE
                    }
            }
        data.aov.s = ifelse( BIG == TRUE, .dobigaov(data$genos,chunksize),
            ifelse(USPARSE==FALSE,summary(aov(lm(trait ~ ., data = data$genos))),
                   summary(aov(SparseM::slm(trait ~ ., data=data$genos),data=data$genos))))
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
                        stop("rsq == NA found")
                    }
                ## adj. r^2 due just to mutations at this freq. bin
                adj.rsq[IDX] =  1 - ( (sum.sq[n] + sum(sum.sq[which(alleleCounts != ac)]))/ttl.sum.sq )*(sum(DF)/(DF[n] + length(which(alleleCounts != ac))))
                if(is.na(adj.rsq[IDX]))
                    {
                        stop("adj.rsq == NA found")
                    }
                IDX=IDX+1
            }
        rsq = cumsum(rsq)
        adj.rsq = cumsum(adj.rsq)
        return(.fillvpv1matrix(data.frame(p,rsq,adj.rsq),twoN))
    }
