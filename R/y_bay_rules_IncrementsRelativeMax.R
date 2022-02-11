# --------------------------------------------------
# Increments control based on relative differences in intervals
# Difference to original function:
# !!using max given dose as reference!!
# Probably the original function can be modifies with an extra parameter
# --------------------------------------------------

#' Increments control based on relative differences in intervals
#'
#' Note that \code{intervals} is to be read as follows. If for example,
#' we want to specify three intervals: First 0 to less than 50, second at least
#' 50 up to less than 100 mg, and third at least 100 mg, then we specify
#' \code{intervals} to be \code{c(0, 50, 100)}. That means, the right
#' bound of the intervals are exclusive to the interval, and the last interval
#' goes from the last value until infinity.
#'
#' @slot intervals a vector with the left bounds of the relevant intervals
#' @slot increments a vector of the same length with the maximum allowable
#' relative increments in the \code{intervals}
#' 
#' @example examples/Rules-class-IncrementsRelative.R
#' @export
#' @keywords classes
.IncrementsRelativeMax <-
  setClass(Class="IncrementsRelativeMax",
           representation(intervals="numeric",
                          increments="numeric"),
           prototype(intervals=c(0, 2),
                     increments=c(2, 1)),
           contains="Increments",
           validity=
             function(object){
               o <- Validate()
               
               o$check(identical(length(object@increments),
                                 length(object@intervals)),
                       "increments must have same length as intervals")
               o$check(! is.unsorted(object@intervals, strictly=TRUE),
                       "intervals has to be sorted and have unique values")
               
               o$result()
             })
validObject(.IncrementsRelativeMax())

#' Initialization function for "IncrementsRelative"
#'
#' @param intervals see \code{\linkS4class{IncrementsRelative}}
#' @param increments see \code{\linkS4class{IncrementsRelative}}
#' @return the \code{\linkS4class{IncrementsRelative}} object
#'
#' @export
#' @keywords methods
IncrementsRelativeMax <- function(intervals,
                                  increments)
{
  .IncrementsRelativeMax(intervals=intervals,
                         increments=increments)
}

# --------------------------------------------------
# The maximum allowable relative increments in intervals method
# !!using max given dose as reference!!
# --------------------------------------------------

#' @describeIn maxDose Determine the maximum possible next dose based on
#' relative increments
#' 
#' @example examples/Rules-method-maxDose-IncrementsRelative.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementsRelativeMax",
                      data="Data"),
          def=
            function(increments, data, ...){
              # determine what was the maximum applied dose
              maxDose <- max(data@x)
              
              # determine in which interval this dose was
              maxInterval <-
                findInterval(x=maxDose,
                             vec=increments@intervals)
              
              # so the maximum next dose is
              ret <-
                (1 + increments@increments[maxInterval]) *
                maxDose
              
              return(ret)
            })
