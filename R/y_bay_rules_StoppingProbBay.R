#' Stopping rule based on the posterior probability of the full model. 
#' The probability that a certain dose is toxic or safe will be calculated
#' and compared to a threshold. The study is stopped in case that the
#' probability is below or above the threshold (toxic or safe).
#'
#' @slot target the target toxicity
#' @slot prob the threshold probability
#' @slot doseeval the dose were the rule is evaluated
#' @slot direction the direction, is the dose toxic or safe
#' @slot dosetested must the dose be already tested
#' @slot label label of stopping rule
#' 
#' @export
#' @keywords classes 
.StoppingProbBay <-
  setClass(Class="StoppingProbBay",
           representation(target="numeric",
                          prob="numeric",
                          doseeval="numeric",
                          direction="character",
                          dosetested="logical",
                          label="character"),
           prototype(target=0.33,
                     prob=0.9,
                     doseeval=10,
                     direction="below",
                     dosetested=TRUE,
                     label="first tox model"),
           contains="Stopping")


#' Initialization function for "StoppingProbBay"
#'
#' @param target see \code{\linkS4class{StoppingProbBay}}
#' @param prob see \code{\linkS4class{StoppingProbBay}}
#' @param doseeval see \code{\linkS4class{StoppingProbBay}}
#' @param direction see \code{\linkS4class{StoppingProbBay}}
#' @param dosetested see \code{\linkS4class{StoppingProbBay}}
#' @param label see \code{\linkS4class{StoppingProbBay}}
#'
#' @export
#' @keywords methods 
StoppingProbBay<- function(target,
                           prob,
                           doseeval,
                           direction,
                           dosetested = TRUE,
                           label)
{
  .StoppingProbBay(target=target,
                   prob=prob,
                   doseeval=doseeval,
                   direction=direction,
                   dosetested=dosetested,
                   label=label)
}


#' @describeIn stopTrial Stopping rule based on posterior probability
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingProbBay",
                      dose="numeric",
                      samples="Samples",
                      model="Model",
                      data="ANY"),
          def=
            function(stopping, dose,
                     samples, model, data, ...){
              # first we have to get samples from the dose-tox
              # curve at the dose.
              probSamples <- prob(dose=stopping@doseeval,
                                  model,
                                  samples)
              
              # what is the probability that the Dose is toxic/safe?
              if      (stopping@direction=='safe') {prob <- mean(probSamples < stopping@target)}
              else if (stopping@direction=='toxic'){prob <- mean(probSamples > stopping@target)}
              
              # Prob high enough?
              probreached <- prob >= stopping@prob 
              
              #dose tested?
              #initialize tested
              tested <- ifelse(stopping@doseeval %in% data@x, TRUE, FALSE)
              
              # so can we stop?
              doStop <- ifelse(stopping@dosetested, probreached & tested, probreached) 
              
              # generate message
              text <-
                paste(doStop, ' : ', stopping@label,
                      " : Probability of dose ",
                      stopping@doseeval,
                      " beeing ",
                      stopping@direction,
                      " is ",
                      round(prob * 100),
                      "%, thus ",
                      ifelse(probreached, "above", "below"),
                      " the required ",
                      round(stopping@prob * 100),
                      "%",
                      if(stopping@dosetested){
                        paste(", and dose is", ifelse(tested, " "," not "), "tested", sep="")
                      }
                      , sep="")
              
              # return both
              return(structure(doStop,
                               message=text))
            })
