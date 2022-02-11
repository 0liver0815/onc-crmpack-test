# --------------------------------------------------
# Increments control based on relative differences in terms of DLTs
# Difference to original function:
# !!using max given dose as reference!!
# Probably the original function can be modifies with an extra parameter
# --------------------------------------------------

#' Increments control based on relative differences in terms of DLTs
#' using max given dose as reference
#'
#' Note that \code{DLTintervals} is to be read as follows. If for example,
#' we want to specify three intervals: First 0 DLTs, second 1 or 2 DLTs, and
#' third at least 3 DLTs, then we specify
#' \code{DLTintervals} to be \code{c(0, 1, 3)}. That means, the right
#' bound of the intervals are exclusive to the interval -- the vector only
#' gives the left bounds of the intervals. The last interval goes from 3 to
#' infinity.
#'
#' @slot DLTintervals an integer vector with the left bounds of the relevant
#' DLT intervals
#' @slot increments a vector of the same length with the maximum allowable
#' relative increments in the \code{DLTintervals}
#'
#' @example examples/Rules-class-IncrementsRelativeDLT.R
#' @export
#' @keywords classes
.IncrementsRelativeDLTmax <-
  setClass(Class="IncrementsRelativeDLTmax",
           representation(DLTintervals="integer",
                          increments="numeric"),
           prototype(DLTintervals=as.integer(c(0, 1)),
                     increments=c(2, 1)),
           contains="Increments",
           validity=
             function(object){
               o <- Validate()
               
               o$check(identical(length(object@increments),
                                 length(object@DLTintervals)),
                       "increments must have same length as DLTintervals")
               o$check(! is.unsorted(object@DLTintervals, strictly=TRUE),
                       "DLTintervals has to be sorted and have unique values")
               o$check(all(object@DLTintervals >= 0),
                       "DLTintervals must only contain non-negative integers")
               
               o$result()
             })
validObject(.IncrementsRelativeDLTmax())


#' Initialization function for "IncrementsRelativeDLTmax"
#'
#' @param DLTintervals see \code{\linkS4class{IncrementsRelativeDLTmax}}
#' @param increments see \code{\linkS4class{IncrementsRelativeDLTmax}}
#' @return the \code{\linkS4class{IncrementsRelativeDLTmax}} object
#'
#' @export
#' @keywords methods
IncrementsRelativeDLTmax <- function(DLTintervals,
                                     increments)
{
  .IncrementsRelativeDLTmax(DLTintervals=safeInteger(DLTintervals),
                            increments=increments)
}


# --------------------------------------------------
# The maximum allowable relative increments in terms of DLTs
# !!using max given dose as reference!!
# --------------------------------------------------

#' @describeIn maxDose Determine the maximum possible next dose based on
#' relative increments determined by DLTs so far
#' 
#' @example examples/Rules-method-maxDose-IncrementsRelativeDLT.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementsRelativeDLTmax",
                      data="Data"),
          def=
            function(increments, data, ...){
              # determine what was the maximum applied dose
              maxDose <- max(data@x)
              
              # determine how many DLTs have occurred so far
              dltHappened <- sum(data@y)
              
              # determine in which interval this is
              interval <-
                findInterval(x=dltHappened,
                             vec=increments@DLTintervals)
              
              # so the maximum next dose is
              ret <-
                (1 + increments@increments[interval]) *
                maxDose
              
              return(ret)
            })