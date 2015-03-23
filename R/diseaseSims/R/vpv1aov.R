.fillvpv1matrix = function(mm,twoN)
    {
        rv = matrix(data=NA,ncol=3,nrow=twoN)
        rv[,1] = (1:twoN)/twoN
        for( i in 1:nrow(mm) )
            {
                P = mm[i,1]*twoN
                rv[P,2]=mm[i,2]
                rv[P,3]=mm[i,3]
            }
        return (rv);
    }

#'Fun with aov
vpv1aov = function(data, ofilename , replicateID, append = FALSE)
    {
        ##Fit the model and summarize
        ##We coerce the matrix to a data frame so that we can get sum of squares, etc.,
        ##per marker
        data.aov.s = summary(aov(SparseM::slm(data$trait ~ ., data=as.data.frame(data$genos)),data=as.data.frame(data$genos)))

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
