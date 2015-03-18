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
        if( modelname != "recessive" |
           modelname != "additive" |
           modelname != "multiplicative" |
           modelname != "popgen" )
            {
                stop("Error: invalid model name.  Must be one of: recessive, additive, multiplicative, or popgen")
            }
        return (.getRiskVariantMatrixDetails(model,popfilename,popfileOffset,recordID,dominance))
    }
