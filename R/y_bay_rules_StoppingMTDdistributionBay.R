# Note: this function is similar to StoppingProbBay
# i.e. it should deliver the same results
#
# This function takes the target tox and calculates the dose value that
# corresponds to the target toxicity. The dose value is compared to the
# specified dose.
#
# The function StoppingProbBay takes a dose level and calculates the toxicity
# probability of the dose level. The toxicicty probability is compared to the
# target toxcicity.

# -----------------------------------------------------------------------------
# Implement a stopping rule based on the posterior distribution
# It can be stopped because a certain dose is toxic or safe
# The probability that a dose is toxic is calculated and the trial is stopped
# if the probability is below or above the threshold 
# -----------------------------------------------------------------------------

#' stopping rule based on the posterior distribution 
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
.StoppingMTDdistributionBay <-
  setClass(Class="StoppingMTDdistributionBay",
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


#' Initialization function for "StoppingMTDdistributionBay"
#'
#' @param target see \code{\linkS4class{StoppingMTDdistributionBay}}
#' @param prob see \code{\linkS4class{StoppingMTDdistributionBay}}
#' @param doseeval see \code{\linkS4class{StoppingMTDdistributionBay}}
#' @param direction see \code{\linkS4class{StoppingMTDdistributionBay}}
#' @param dosetested see \code{\linkS4class{StoppingMTDdistributionBay}}
#' @param label see \code{\linkS4class{StoppingMTDdistributionBay}}
#'
#' @export
#' @keywords methods 
StoppingMTDdistributionBay<- function(target,
                                      prob,
                                      doseeval,
                                      direction,
                                      dosetested = TRUE,
                                      label)
{
  .StoppingMTDdistributionBay(target=target,
                              prob=prob,
                              doseeval=doseeval,
                              direction=direction,
                              dosetested=dosetested,
                              label=label)
}


#' @describeIn stopTrial Stopping rule based on posterior dist
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingMTDdistributionBay",
                      dose="numeric",
                      samples="Samples",
                      model="Model",
                      data="ANY"),
          def=
            function(stopping, dose,
                     samples, model, data, ...){
              # First, generate the MTD samples.
              
              # add prior data and samples to the
              # function environment so that they
              # can be used.
              mtdSamples <- dose(prob=stopping@target,
                                 model,
                                 samples)
              
              # what is the probability that the MTD is below/above this dose?
              if      (stopping@direction=='below'){prob <- mean(mtdSamples < stopping@doseeval)}
              else if (stopping@direction=='above'){prob <- mean(mtdSamples > stopping@doseeval)}
              
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
                      " : Probability of MTD ",
                      stopping@direction,
                      " dose ",
                      stopping@doseeval,
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

