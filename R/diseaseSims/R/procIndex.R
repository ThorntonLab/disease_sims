#' Calculate offset values from index files
#' @param index A data frame or matrix representing the index file data.
#' @param column  The relevant column in index.
#' @param recordno The record number whose index we want to look up
#' @details
#' The programs in disease_sims use file locking to allow multiple instances
#' of a program to write output to the same file.  In order to keep track
#' of particular replicates across various output files, an index file is needed.
#' These index files consist of at least two columns.  The first column is
#' an integer identifying a specific replicate.  (Typically, this number is the
#' same as the task ID number used in the array job to generate the data.)
#' The remaining columns relate to either the offset (in bytes) where the relevant
#' record starts (in the case of a file generated without gzip compression) or
#' the size of the record in bytes (for output generating using gzip compression).
#' This function goes through such indexes and returns the offset for the record,
#' allowing a seek/gzseek to be performed on an input file.
procIndex = function(index,column,recordno)
    {
        if (! is.data.frame(index) && ! is.matrix(index) )
            {
                stop("procIndex error: data must be matrix or data frame with at least two columns")
            }
        if ( ncol(index) < 2 )
            {                
                stop("procIndex error: data frame has fewer than two columns")
            }
        if( column <= 1  || column > ncol(index) )
            {
                stop(paste("procIndex error: column value of",
                           column,
                           "is out of bounds.  Data have",
                           ncol(index),"columns"))
            }
        record = which(index[,1] == recordno)
        if (length(record) == 0)
            {
                stop(paste("procIndex error: record number",
                           recordno,
                           "is not present in the first column of the data"))
            }
        if( record == 1 )
            {
                return(0)
            }
        #check if index values are a monotonically-increasing series.
        #If so, then they are offsets in bytes.
        #If not, then they are sizeof(record) in bytes
        ispos = if (length(which(diff(index[,column])>0)) == nrow(index)-1) TRUE else FALSE
        if( ispos )
            {
                return (index[record,column])
            }
        else
            {
                return( sum(index[1:(record-1),column]) )
            }
    }
