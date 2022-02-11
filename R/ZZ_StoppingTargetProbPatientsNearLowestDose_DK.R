## -----------------------------------------------------------------------------------------
## New addition for [Rules-class.R] 
## -----------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------
## Stopping based on probability of target tox interval and number of patients near to the lowest dose
## ---------------------------------------------------------------------------------------------------

##' Stop based on number of patients near to next best dose
##'
##' @slot target the target toxicity interval, e.g. \code{c(0.2, 0.35)}
##' @slot prob required target toxicity probability (e.g. \code{0.4})
##' for reaching sufficient precision
##' @slot nPatients number of required patients
##' @slot percentage percentage (between 0 and 100) within the next best dose
##' the patients must lie
##' 
##' @keywords classes
##' @export
.StoppingTargetProbPatientsNearLowestDose <-
  setClass(Class="StoppingTargetProbPatientsNearLowestDose",
           representation(target="numeric",
                          prob="numeric",
                          nPatients="integer",
                          percentage="numeric"),
           prototype(target=c(0.2, 0.35),
                     prob=0.4,
                     nPatients=10L,
                     percentage=50),
           contains="Stopping",
           validity=function(object){
             o <- Validate()
             
             o$check(is.probRange(object@target),
                     "target must be probability range")
             o$check(is.probability(object@prob,
                                    bounds=FALSE),
                     "prob must be probability > 0 and < 1")
             o$check((object@nPatients > 0L) && is.scalar(object@nPatients),
                     "nPatients must be positive scalar")
             o$check(is.probability(object@percentage / 100),
                     "percentage must be between 0 and 100")
             
             o$result()
           })
validObject(.StoppingTargetProbPatientsNearLowestDose())


##' Initialization function for "StoppingTargetProbPatientsNearLowestDose"
##'
##' @param target see \code{\linkS4class{StoppingTargetProb}}
##' @param prob see \code{\linkS4class{StoppingTargetProb}}
##' @return the \code{\linkS4class{StoppingTargetProb}} object
##' @param nPatients see \code{\linkS4class{StoppingPatientsNearDose}}
##' @param percentage see \code{\linkS4class{StoppingPatientsNearDose}}
##' @return the \code{\linkS4class{StoppingPatientsNearDose}} object
##'
##' @export
##' @keywords methods
StoppingTargetProbPatientsNearLowestDose <- function(target,
                                                     prob,
                                                     nPatients,
                                                     percentage)
{
  .StoppingTargetProbPatientsNearLowestDose(target=target,
                                            prob=prob,
                                            nPatients=safeInteger(nPatients),
                                            percentage=percentage)
}

## -----------------------------------------------------------------------------------------
## New addition for [Rules-methods.R] 
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------------
## Stopping based on probability of target tox interval and number of patients near to the lowest dose
## -----------------------------------------------------------------------------------------------

##' @describeIn stopTrial Stop based on probability of target tox interval and 
##' number of patients near to the lowest dose

setMethod("stopTrial",
          signature=
            signature(stopping="StoppingTargetProbPatientsNearLowestDose",
                      dose="numeric",
                      samples="Samples",
                      model="Model",
                      data="ANY"),
          def=
            function(stopping, dose, samples, model, data, ...){
              ## first we have to get samples from the dose-tox
              ## curve at the dose.
              dose <- data@doseGrid[1]
              probSamples <- prob(dose=dose,
                                  model,
                                  samples)
              
              ## Now compute probability to be in target interval
              probTarget <-
                mean((probSamples >= stopping@target[1]) &
                       (probSamples <= stopping@target[2]))
              
              ## so can we stop due to toxicity interval?
              doStopprobTarget <- probTarget >= stopping@prob
              
              
              lower <- (100 - stopping@percentage) / 100 * dose
              upper <- (100 + stopping@percentage) / 100 * dose
              
              ## how many patients lie there?
              nPatients <- sum((data@x >= lower) & (data@x <= upper))
              
              ## so can we stop due to sufficient number of patients?
              doStopnPatients <- nPatients >= stopping@nPatients
              
              ## so can we stop overall?
              doStop <- all(c(doStopprobTarget,doStopnPatients)=="TRUE")
              
              ## generate message
              text <-
                paste("Probability for target toxicity is",
                      round(probTarget * 100),
                      " % for lowest dose ",
                      dose,
                      " and thus ",
                      ifelse(doStopprobTarget, "above", "below"),
                      " the required ",
                      round(stopping@prob * 100),
                      "%"," and ",
                      nPatients,
                      " patients lie within ",
                      stopping@percentage,
                      " % of the lowest dose ",
                      dose,
                      ". This ",
                      ifelse(doStopnPatients, "reached", "is below"),
                      " the required ",
                      stopping@nPatients,
                      " patients",
                      sep="")
              
              ## return both
              return(structure(doStop,
                               message=text))
            })

