.fillvpv1matrix = function(mm,twoN)
    {
        rv = data.frame((1:twoN)/twoN,NA,NA,NA,NA)
        colnames(rv) <- c("V1","V2","V3","V4","V5")
        for( i in 1:nrow(mm) )
            {
                P = mm[i,1]*twoN
                rv[P,2]=mm[i,2]
                rv[P,3]=mm[i,3]
                rv[P,4]=mm[i,4]
                rv[P,5]=mm[i,5]
            }
        ##Courtesy of http://stackoverflow.com/questions/23340150/using-dplyr-window-functions-to-make-trailing-values
        rv %>%
            mutate(dummy = cumsum(0 + !is.na(V2))) %>%
                group_by(dummy,add=TRUE) %>%
                    mutate(filled = nth(V2,1)) %>%
                        mutate(filled2 = nth(V3,1)) %>%
                            mutate(filled3 = nth(V4,1)) %>%
                                mutate(filled4 = nth(V5,1)) %>%
                        ungroup() %>%
                            select(-V2,-V3,-V4,-V5,-dummy) -> m2
        colnames(m2) <- c("p","rsq.A","adj.rsq.A","rsq.D","adj.rsq.D")
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
                BIGGIE <- update( BIGGIE, moredata=data[IIa:IIb,] )
                IIa = IIb + 1
                IIb = min(IIa + chunksize - 1, FINAL)
            }
        RV =summary(aov(BIGGIE,data=data))
        return(RV)
    }

#' Fit linear models to genotypes
#' @param data A list returned by getRiskVariantMatrix
#' @param method Method to use for model fitting.  Must be one of "biglm","slm", or "lm", corresponding to the biglm, SparseM, and base packages, respectively.
#' @param chunksize If the data set is massive, we'll use biglm::biglm for the regression, breaking data up into sets of chunksize rows.  In our testing, larger chunk sizes result in lower run-times.
#' @return A data frame with 3 columns: allele frequency, and then estimates of the total variance in the trait explained by markers are <= that frequency.  The estimates make up the last two columns, and are based on the R^2 and adjusted R^2 of the linear model, respectively.
vpv1aovDomi = function(data, method = "lm", chunksize=5000)
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
        ROWS.A = names(data$genos[seq(2,ncol(data$genos)-1,2)])
        ROWS.D = names(data$genos[seq(3,ncol(data$genos),2)])

        ROWS =sapply(as.array(rownames(data.aov.s[[1]])[1:(length(rownames(data.aov.s[[1]]))-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)
        
        ROWS.R.A= ROWS[ROWS%in% ROWS.A]

        ROWS.R.D=ROWS[ROWS%in% ROWS.D]
        
        ##DOMINANCE MATRIX DOESN'T PLAY NICE FOR ALLELE FREQUENICES
        mutc = colSums(data$genos[,2:ncol(data$genos)])
        #set the count of the dominance component to the count of the corresponding additive component
        mutc[seq(2,length(mutc),2)] = mutc[seq(1,length(mutc)-1,2)]
        
        alleleCounts.A = as.integer(mutc[ROWS.R.A])
        alleleCounts.D = as.integer(mutc[ROWS.R.D])
        A.index = which( ROWS %in% ROWS.R.A )
        D.index = which( ROWS %in% ROWS.R.D )
           #alleleCounts= colSums(data$genos[,sapply(as.array(ROWS[1:(length(ROWS)-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)]))

        ##Initalize the matrix to return.
        ##The columns will be: mutant allele frequency, r^2, adj. r^2
      
        
        sum.sq = data.aov.s[[1]]$'Sum Sq'
        ttl.sum.sq = sum(sum.sq)
        DF = data.aov.s[[1]]$'Df'
        n = length(sum.sq)
        RSS = sum.sq[n]
        ##Populate the matrix, starting with the rares
        for(component in c("A","D")){
            IDX=1
            alleleCounts =if(component=="A"){alleleCounts.A}else{alleleCounts.D}
            index =if(component=="A"){A.index}else{D.index}
            index.not= if(component=="A"){D.index}else{A.index}
            p.temp = array( dim = length(unique(alleleCounts)) )
            rsq.temp = array( dim = length(unique(alleleCounts)) )
            adj.rsq.temp = array( dim = length(unique(alleleCounts)) )
            SS.temp = sum.sq[index]
            SS.not.temp = sum.sq[index.not]
            for( ac in sort(unique(alleleCounts)) )
                {
                    ac.sum.sq = SS.temp[which(alleleCounts == ac)]
                    p.temp[IDX] = ac/twoN
                    ## r^2 due just to mutation at this freq. bin
                    rsq.temp[IDX] = sum(ac.sum.sq)/ttl.sum.sq
                    if(is.na(rsq.temp[IDX]))
                        {
                        stop("rsq == NA found")
                    }
                    ## adj. r^2 due just to mutations at this freq. bin
                    adj.rsq.temp[IDX] =  1 - ( (RSS + sum(SS.temp[which(alleleCounts != ac)]) + sum(SS.not.temp) )/ttl.sum.sq)*(sum(DF)/(DF[n] + length(which(alleleCounts != ac)) + length(SS.not.temp)))
                    if(is.na(adj.rsq.temp[IDX]))
                        {
                            stop("adj.rsq == NA found")
                        }
                    IDX=IDX+1
                }
            if( component == "A")
                {
                    p.A =  p.temp
                    rsq.A =  cumsum(rsq.temp)
                    adj.rsq.A =  cumsum(adj.rsq.temp)
                }
            else if (component =="D")
                {
                    p.D = c(0, p.temp)
                    rsq.D =  c(0,cumsum(rsq.temp))
                    adj.rsq.D =  c(0,cumsum(adj.rsq.temp))
                }
            
        }

       
        rsq.D.A =  rsq.D[sapply(p.A,function(x)max(which(p.D<=x)))]
        adj.rsq.D.A =  adj.rsq.D[sapply(p.A,function(x)max(which(p.D<=x)))]
     
        return(.fillvpv1matrix(data.frame(p.A,rsq.A,adj.rsq.A,rsq.D.A,adj.rsq.D.A),twoN))
    }
