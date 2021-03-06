#' Return a genotype matrix based on risk variants with dominance from a population subset
#' @param popfilename The (binary, gzipped) file containing a population simulated by TFL2013_ind
#' @param popfileOffset The offset (in bytes) of the relevant record in popfilename
#' @param modelName The model of gene action.  Must be one of: recessive, additive, multiplicative, or popgen
#' @param dominance The dominance of a risk mutation.  Only applies to the "popgen" model of gene action
#' @param selectedOnly True or False, get causative sites only
#' @param subset A vector of 0/1, specifying which individuals to take from the population
#' @param nsample An integer value specifying the size of the subset being taken
#' @return genos, which is an data frame of genotypes. Rows = individuals. The first column contains trait
#' values for each individual, and the remaining columns = 0,1,2 copies of risk mutation
#' @return esizes, which is a numeric vector of effect sizes for each mutation in genos
#' @return positions, which are the positions of every matrix in genos
#' @details
#' The order of the columns in "genos" is in descending order of both frequency and absolute value of effect size
#' If phenofilename and phenofileOffset are both provided, then "traits" corresponds to individual phenotypes.
#' If these two arguments are NOT provided, then "traits" corresponds to the genetic component of trait value,
#' which is calculated by the parameter based to modelName.
#' Risk variant frequencies are calculated by colSums(genos)[-1]/(2*nrow(genos)), in case you want to filter
getVariantMatrixDominanceSubset = function(popfilename,
    popfileOffset,
    modelName = "recessive",
    dominance = 0,
    selectedOnly = TRUE,
    subset,
    nsample)
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
        if(modelName != "recessive" &
           modelName != "additive" &
           modelName != "multiplicative" &
           modelName != "popgen")
            {
                stop("Error: invalid model name.  Must be one of: recessive, additive, multiplicative, or popgen")
            }
        
        XX = .getVariantMatrixDominanceDetailsSubset(modelName,popfilename,popfileOffset,subset,nsample,dominance,selectedOnly)
        return (XX);
    }

#' Return a genotype matrix based on risk variants from a population
#' @param popfilename The (binary, gzipped) file containing a population simulated by TFL2013_ind
#' @param popfileOffset The offset (in bytes) of the relevant record in popfilename
#' @param modelName The model of gene action.  Must be one of: recessive, additive, multiplicative, or popgen
#' @param dominance The dominance of a risk mutation.  Only applies to the "popgen" model of gene action
#' @param subset A vector of 0/1, specifying which individuals to take from the population
#' @param nsample An integer value specifying the size of the subset being taken
#' @return genos, which is an data frame of genotypes. Rows = individuals. The first column contains trait
#' values for each individual, and the remaining columns = 0,1,2 copies of risk mutation
#' @return esizes, which is a numeric vector of effect sizes for each mutation in genos
#' @return positions, which are the positions of every matrix in genos
#' @details
#' The order of the columns in "genos" is in descending order of both frequency and absolute value of effect size
#' If phenofilename and phenofileOffset are both provided, then "traits" corresponds to individual phenotypes.
#' If these two arguments are NOT provided, then "traits" corresponds to the genetic component of trait value,
#' which is calculated by the parameter based to modelName.
#' Risk variant frequencies are calculated by colSums(genos)[-1]/(2*nrow(genos)), in case you want to filter
#'
#' This function simply calls getVariantMatrix with selectedOnly = TRUE
getRiskVariantMatrixDominanceSubset = function(popfilename,
    popfileOffset,
    modelName = "recessive",
    dominance = 0,
    subset,
    nsample)
    {
    return( getVariantMatrixDominanceSubset(popfilename,popfileOffset,modelName,dominance,TRUE,subset,nsample) )
}
