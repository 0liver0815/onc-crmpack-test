## --------------------------------------------------
## Stopping based on number of patients near to next best dose
## No change to original rule other than output text
## --------------------------------------------------

##' Stop based on number of patients near to next best dose
##' !! No rule change to standard StoppingPatientsNearDose
##' Only change is returned text, when h stopping rule is triggered
##'
##' @slot nPatients number of required patients
##' @slot percentage percentage (between 0 and 100) within the next best dose
##' the patients must lie
##' @slot label label of stopping rule
##' 
##' @example examples/Rules-class-StoppingPatientsNearDose.R
##' @keywords classes
##' @export
.StoppingPatientsNearDose2 <-
  setClass(Class="StoppingPatientsNearDose2",
           representation(label="character"),
           prototype(label="N9"),
           contains="StoppingPatientsNearDose")

##' Initialization function for "StoppingPatientsNearDose2"
##'
##' @param nPatients see \code{\linkS4class{StoppingPatientsNearDose}}
##' @param percentage see \code{\linkS4class{StoppingPatientsNearDose}}
##' @param label see \code{\linkS4class{StoppingPatientsNearDose}}
##' @return the \code{\linkS4class{StoppingPatientsNearDose}} object
##'
##' @export
##' @keywords methods
StoppingPatientsNearDose2 <- function(nPatients,
                                      percentage,
                                      label)
{
  .StoppingPatientsNearDose2(nPatients=safeInteger(nPatients),
                             percentage=percentage,
                             label=label)
} 

## -------------------------------------------------------------
## Stopping based on number of patients near to next best dose
## No change to original rule other than returned text
## -------------------------------------------------------------

##' @describeIn stopTrial Stop based on number of patients near to next best
##' dose. Identical behavior to the original rule StoppingPatientsNearDose
##' The only change is that additional text is putted, so that it can be
##' identified if the rule is triggered
##' 
##' @example examples/Rules-method-stopTrial-StoppingPatientsNearDose.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingPatientsNearDose2",
                      dose="numeric",
                      samples="ANY",
                      model="ANY",
                      data="Data"),
          def=
            function(stopping, dose, samples, model, data, ...){
              ## determine the range where the cohorts must lie in
              lower <- (100 - stopping@percentage) / 100 * dose
              upper <- (100 + stopping@percentage) / 100 * dose
              
              ## how many patients lie there?
              nPatients <- sum((data@x >= lower) & (data@x <= upper))
              
              ## so can we stop?
              doStop <- nPatients >= stopping@nPatients
              
              ## generate message
              text <- paste(doStop, ' : ', stopping@label,
                            ' : ', nPatients,
                            " patients lie within ",
                            stopping@percentage,
                            "% of the next best dose ",
                            dose,
                            ". This ",
                            ifelse(doStop, "reached", "is below"),
                            " the required ",
                            stopping@nPatients,
                            " patients",
                            sep="")
              
              ## return both
              return(structure(doStop,
                               message=text))
            }) 
