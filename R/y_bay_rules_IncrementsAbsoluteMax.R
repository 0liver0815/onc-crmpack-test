# --------------------------------------------------
# Increments control based on absolute differences in intervals
# Difference to original function:
# !!using max given dose as reference!!
# Probably the original function can be modifies with an extra parameter
# --------------------------------------------------

#' Increments control based on absolute differences in intervals
#'
#' Note that \code{intervals} is to be read as follows. If for example,
#' we want to specify three intervals: First 0 to less than 50, second at least
#' 50 up to less than 100 mg, and third at least 100 mg, then we specify
#' \code{intervals} to be \code{c(0, 50, 100)}. That means, the right
#' bound of the intervals are exclusive to the interval, and the last interval
#' goes from the last value until infinity.
#'
#' @slot intervals a vector with the left bounds of the relevant intervals
#' @slot absoluteIncrements a vector of the same length with the maximum allowable
#' relative increments in the \code{intervals}
#' 
#' @example examples/Rules-class-IncrementsRelative.R
#' @export
#' @keywords classes
.IncrementsAbsoluteMax <-
  setClass(Class="IncrementsAbsoluteMax",
           representation(intervals="numeric",
                          absoluteIncrements="numeric"),
           prototype(intervals=c(0, 2),
                     absoluteIncrements=c(0.5, 1)),
           contains="Increments",
           validity=
             function(object){
               o <- Validate()
               
               o$check(identical(length(object@absoluteIncrements),
                                 length(object@intervals)),
                       "increments must have same length as intervals")
               o$check(! is.unsorted(object@intervals, strictly=TRUE),
                       "intervals has to be sorted and have unique values")
               
               o$result()
             })
validObject(.IncrementsAbsoluteMax())

#' Initialization function for "IncrementsAbsoluteMax"
#'
#' @param intervals see \code{\linkS4class{IncrementsAbsoluteMax}}
#' @param absoluteIncrements see \code{\linkS4class{IncrementsAbsoluteMax}}
#' @return the \code{\linkS4class{IncrementsAbsoluteMax}} object
#'
#' @export
#' @keywords methods
IncrementsAbsoluteMax <- function(intervals,
                                  absoluteIncrements)
{
  .IncrementsAbsoluteMax(intervals=intervals,
                         absoluteIncrements=absoluteIncrements)
}

# --------------------------------------------------
# The maximum allowable relative increments in intervals method
# !!using max given dose as reference!!
# --------------------------------------------------

#' @describeIn maxDose Determine the maximum possible next dose based on
#' absolute increments
#' 
#' @example examples/Rules-method-maxDose-IncrementsRelative.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementsAbsoluteMax",
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
                increments@absoluteIncrements[maxInterval] +
                maxDose
              
              return(ret)
            })
