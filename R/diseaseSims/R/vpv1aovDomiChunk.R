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

.dobigaovChunk = function(popfile, offset, model, dominance,popsize,chunksize = 5000)
    {
        FINAL = as.numeric(popsize)
        IIa=1
        IIb=min(chunksize,FINAL)
        ##Make the model expression
        subset = rep(0,popsize)
        subset[IIa:IIb]=1
        print("about to read data")
        data = getRiskVariantMatrixDominanceSubset(popfile,offset,model,dominance,subset,sum(subset))$genos

        #we have to get the allele counts as we read the data to avoid reading it back in later
        mutc = colSums(data[,2:ncol(data)])
        #set the count of the dominance component to the count of the corresponding additive component
       mutc[seq(2,length(mutc),2)] = mutc[seq(1,length(mutc)-1,2)]
        
        alleleCounts = as.integer(mutc)
       
        FORMULA <- paste(colnames(data)[1], "~", paste(colnames(data)[-1], collapse=" + "))
        print("about to fit big lm")
        BIGGIE = biglm::biglm( as.formula(FORMULA) ,data = data[IIa:IIb,] )
        subset[IIa:IIb]=0
        IIa = IIb + 1
        IIb = min(IIa + chunksize - 1, FINAL)
        while( as.numeric(IIa) < FINAL )
            {
               
                subset[IIa:IIb]=1
                print(IIa)
                print(IIb)
       

                data = getRiskVariantMatrixDominanceSubset(popfile,offset,model,dominance,subset,sum(subset))$genos

                mutc = colSums(data[,2:ncol(data)])
                mutc[seq(2,length(mutc),2)] = mutc[seq(1,length(mutc)-1,2)]
                alleleCounts = alleleCounts + as.integer(mutc)

                BIGGIE =  update( BIGGIE, moredata=data)
                subset[IIa:IIb]=0
                IIa = IIb + 1
                IIb = min(IIa + chunksize - 1, FINAL)
            }
        #data = getRiskVariantMatrixDominance(popfile,offset,model,dominance)$genos
        #SA =summary(aov(BIGGIE,data=data))
        RV =list(SS=c(BIGGIE$qr$D*BIGGIE$qr$thetab^2,BIGGIE$qr$ss)[-1],ac=alleleCounts)

        return(RV)
    }

#' Fit linear models to genotypes
#' @param popfile population file
#' @param offset population file offset
#' @param model name of genetic model
#' @param dominance degree of dominance
#' @param method Method to use for model fitting.  Must be one of "biglm","slm", or "lm", corresponding to the biglm, SparseM, and base packages, respectively.
#' @param chunksize If the data set is massive, we'll use biglm::biglm for the regression, breaking data up into sets of chunksize rows.  In our testing, larger chunk sizes result in lower run-times.
#' @return A data frame with 3 columns: allele frequency, and then estimates of the total variance in the trait explained by markers are <= that frequency.  The estimates make up the last two columns, and are based on the R^2 and adjusted R^2 of the linear model, respectively.
vpv1aovDomiChunk = function(popfile, offset, model, dominance, popsize, method = "biglm", chunksize=5000)
    {
      
        if ( requireNamespace("biglm",quietly=TRUE) == FALSE )
            {
                stop("Error:biglm is required")
            }
        
        RV =  .dobigaovChunk(popfile, offset, model, dominance, popsize, chunksize)
        SS = RV$SS
       
        twoN = 2*(popsize)
        alleleCounts = RV$ac
        alleleCounts.A = alleleCounts[seq(1,length(alleleCounts),2)]
        alleleCounts.D = alleleCounts[seq(2,length(alleleCounts),2)]
        
        ##Get the counts of risk allele frequencies at each marker
        ##GIVEN THAT THE MARKER WERE USED IN THE REGRESSION
        cA = seq(1,length(SS)-1,2)
        cD = seq(2,length(SS),2)

        #ROWS =sapply(as.array(rownames(data.aov.s[[1]])[1:(length(rownames(data.aov.s[[1]]))-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)
        
        #ROWS.R.A= ROWS[ROWS%in% ROWS.A]

        #ROWS.R.D=ROWS[ROWS%in% ROWS.D]
        
        ##DOMINANCE MATRIX DOESN'T PLAY NICE FOR ALLELE FREQUENICES

        #A.index = which( ROWS %in% ROWS.R.A )
        #D.index = which( ROWS %in% ROWS.R.D )
           #alleleCounts= colSums(data$genos[,sapply(as.array(ROWS[1:(length(ROWS)-1)]),function(x) gsub(" ","",x),USE.NAMES=FALSE)]))

        ##Initalize the matrix to return.
        ##The columns will be: mutant allele frequency, r^2, adj. r^2
      
        
        #sum.sq = data.aov.s[[1]]$'Sum Sq'
        TSS = sum(SS)
        n = length(SS)
        RSS = SS[n]
        ##Populate the matrix, starting with the rares
        for(component in c("A","D")){
            IDX=1
            alleleCounts =if(component=="A"){alleleCounts.A}else{alleleCounts.D}
            index =if(component=="A"){cA}else{cD}
            index.not= if(component=="A"){cD}else{cA}
            p.temp = array( dim = length(unique(alleleCounts)) )
            rsq.temp = array( dim = length(unique(alleleCounts)) )
            adj.rsq.temp = array( dim = length(unique(alleleCounts)) )
            SS.temp = SS[index]
            SS.not.temp =SS[index.not]
            for( ac in sort(unique(alleleCounts)) )
                {
                    ac.sum.sq = SS.temp[which(alleleCounts == ac)]
                    p.temp[IDX] = ac/twoN
                    ## r^2 due just to mutation at this freq. bin
                    rsq.temp[IDX] = sum(ac.sum.sq)/TSS
                    if(is.na(rsq.temp[IDX]))
                        {
                        stop("rsq == NA found")
                    }
                    ## adj. r^2 due just to mutations at this freq. bin
                    adj.rsq.temp[IDX] =   1 - ( (RSS + sum(SS.temp[which(alleleCounts != ac)]) + sum(SS.not.temp) )/TSS)*((popsize-1)/(popsize - n - 2 + length(which(alleleCounts != ac)) + length(SS.not.temp)))
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
