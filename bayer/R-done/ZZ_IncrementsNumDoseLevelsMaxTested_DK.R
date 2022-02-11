## -----------------------------------------------------------------------------------------
## New addition for [Rules-class.R] 
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## Increments control based on number of dose levels relative to the max dose tested so far 
## -----------------------------------------------------------------------------------------

##' Increments control based on number of dose levels relative to the max dose tested so far
##'
##' @slot maxLevels scalar positive integer for the number of maximum 
##' dose levels to increment for the next dose. It defaults to 1, 
##' which means that no dose skipping is allowed - the next dose 
##' can be maximum one level higher than the max dose tested so far.
##' 
##' @export
##' @keywords classes
.IncrementsNumDoseLevelsMaxTested <-
  setClass(Class="IncrementsNumDoseLevelsMaxTested",
           representation(maxLevels="integer"),
           prototype(maxLevels=1L),
           contains="Increments",
           validity=
             function(object){
               o <- Validate()
               
               o$check(is.scalar(object@maxLevels) && 
                         is.integer(object@maxLevels) && 
                         object@maxLevels > 0,
                       "maxLevels must be scalar positive integer")
               
               o$result()
             })
validObject(.IncrementsNumDoseLevelsMaxTested())

##' Initialization function for "IncrementsNumDoseLevelsMaxTested"
##'
##' @param maxLevels see \code{\linkS4class{IncrementsNumDoseLevelsMaxTested}}
##' @return the \code{\linkS4class{IncrementsNumDoseLevelsMaxTested}} object
##'
##' @export
##' @keywords methods
IncrementsNumDoseLevelsMaxTested <- function(maxLevels=1)
{
  .IncrementsNumDoseLevelsMaxTested(maxLevels=safeInteger(maxLevels))
}

## -----------------------------------------------------------------------------------------
## New addition for [Rules-methods.R] 
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## The maximum allowable number of dose levels method relative to the max dose tested so far
## -----------------------------------------------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' maximum dose levels to increment for the next dose
##' 
setMethod("maxDose",
          signature=
            signature(increments="IncrementsNumDoseLevelsMaxTested",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the level of the max dose tested so far
              maxDoseLevel <- max(data@xLevel)
              
              ## determine the maximum next dose level
              maxNextDoseLevel <- min(length(data@doseGrid),
                                      maxDoseLevel + increments@maxLevels)
              
              ## so the maximum next dose is
              ret <- data@doseGrid[maxNextDoseLevel]
              
              return(ret)
            })
