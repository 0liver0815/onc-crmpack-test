## --------------------------------------------------
## Increments control based on number of dose levels
## !!using max given dose as reference!!
## --------------------------------------------------

#' Increments control based on number of dose levels.
#' Reference is the maximum applied dose.
#'
#' @slot maxLevels scalar positive integer for the number of maximum 
#' dose levels to increment for the next dose. It defaults to 1, 
#' which means that no dose skipping is allowed - the next dose 
#' can be maximum one level higher than maximum given dose.
#' 
#' @example examples_new/Rules-class-IncrementsNumDoseLevelsMax.R
#' @export
#' @keywords classes
.IncrementsNumDoseLevelsMax <-
  setClass(Class="IncrementsNumDoseLevelsMax",
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
validObject(.IncrementsNumDoseLevelsMax())

#' Initialization function for "IncrementsNumDoseLevelsMax"
#'
#' @param maxLevels see \code{\linkS4class{IncrementsNumDoseLevelsMax}}
#' @return the \code{\linkS4class{IncrementsNumDoseLevelsMax}} object
#'
#' @export
#' @keywords methods
IncrementsNumDoseLevelsMax <- function(maxLevels=1)
{
  .IncrementsNumDoseLevelsMax(maxLevels=safeInteger(maxLevels))
}

## --------------------------------------------------
## The maximum allowable number of dose levels method
## !!using max given dose as reference!!
## --------------------------------------------------

#' @describeIn maxDose Determine the maximum possible next dose based on
#' maximum dose levels to increment for the next dose
#' 
#' @example examples_new/Rules-method-maxDose-IncrementsNumDoseLevelsMax.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementsNumDoseLevelsMax",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the maximum applied dose
              maxDoseLevel <- max(data@xLevel)
              
              ## determine the maximum next dose level
              maxNextDoseLevel <- min(data@nGrid,
                                      maxDoseLevel + increments@maxLevels)
              
              ## so the maximum next dose is
              ret <- data@doseGrid[maxNextDoseLevel]
              
              return(ret)
            })

