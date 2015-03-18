#' Return a genotype matrix based on risk variants from a population
getRiskVariantMatrix = function(popfilename,
    popfileOffset,
    recordID = 0,
    modelName = "recessive",
    dominance = 0)
    {
        if( missing(popfilename) ) {
            stop("Error: popfilename required")
        }
        if( missing(popfileOffset) ) {
            stop("Error: popfileOffset required")
        }
        if( popfileOffset < 0 ) {
            stop("Error: popfileOffset must be >= 0")
        }
        if( modelName != "recessive" &
           modelName != "additive" &
           modelName != "multiplicative" &
           modelName != "popgen" )
            {
                stop("Error: invalid model name.  Must be one of: recessive, additive, multiplicative, or popgen")
            }
        return (.getRiskVariantMatrixDetails(modelName,popfilename,popfileOffset,recordID,dominance))
    }
