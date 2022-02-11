
## -----------------------------------------------------------------------------------------
## Add pbeta in importFrom stats in the programm crmPack-package.r
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## New addition for [Rules-class.R] 
## -----------------------------------------------------------------------------------------

## --------------------------------------------------
## Stopping based on the lowest dose meeting the hard safety criteria using Beta based DLT probability
## --------------------------------------------------

##' Stop based on the lowest dose meeting the hard safety criteria using 
##' Beta based DLT probability, i.e. 1- pbeta(target, x+a, n-x+b) < prob
##'
##' @slot target toxicity target for a dose
##' @slot prob probability of dose being toxic
##' @slot a shape parameter a>0 of probability distribution Beta (a,b)
##' @slot b shape parameter b>0 of probability distribution Beta (a,b)
##' 
##' @keywords classes
##' @export

.StoppingLowestDoseHSRBeta_DK <-
  setClass(Class="StoppingLowestDoseHSRBeta_DK",
           representation(target="numeric",
                          prob="numeric",
                          a="numeric",
                          b="numeric"),
           prototype(target=0.3,
                     prob=0.9,
                     a=1,
                     b=1),
           contains="Stopping",
           validity=
             function(object){
               o <- Validate()
               
               o$check(is.probability(object@target,
                                      bounds=FALSE),
                       "target must be probability > 0 and < 1")
               o$check(is.probability(object@prob,
                                      bounds=FALSE),
                       "prob must be probability > 0 and < 1")
               o$check(is.numeric(object@a) & object@a > 0,
                       "Beta distribution shape parameter a must me a real number > 0")
               o$check(is.numeric(object@b) & object@b > 0,
                       "Beta distribution shape parameter b must me a real number > 0")
               o$result()
             })
validObject(.StoppingLowestDoseHSRBeta_DK())

##' Initialization function for "StoppingLowestDoseHSRBeta_DK"
##'
##' @param target see \code{\linkS4class{StoppingLowestDoseHSRBeta_DK}}
##' @param prob see \code{\linkS4class{StoppingLowestDoseHSRBeta_DK}}
##' @param a see \code{\linkS4class{StoppingLowestDoseHSRBeta_DK}}
##' @param b see \code{\linkS4class{StoppingLowestDoseHSRBeta_DK}}
##' @return the \code{\linkS4class{StoppingLowestDoseHSRBeta_DK}} object
##'
##' @export
##' @keywords methods

StoppingLowestDoseHSRBeta_DK <- function(target,prob,a,b)
{
  .StoppingLowestDoseHSRBeta_DK(target=target,
                             prob=prob,
                             a=a,
                             b=b)
}

## -----------------------------------------------------------------------------------------
## New addition for [Rules-methods.R] 
## -----------------------------------------------------------------------------------------

## -------------------------------------------------------------
## Stopping based on the lowest dose meeting the hard safety criteria using Beta based DLT probability
## -------------------------------------------------------------

##' @describeIn stopTrial Stop based on the lowest dose meeting the hard 
##' safety criteria using Beta based DLT probability

setMethod("stopTrial",
          signature=
            signature(stopping="StoppingLowestDoseHSRBeta_DK",
                      dose="numeric",
                      samples="ANY",
                      model="ANY",
                      data="Data"),
          def=
            function(stopping, dose, samples, model, data, ...){
              #Load data dose vector and find the number of total subjects dosed at this the lowest planned dose
              dose <- data@doseGrid[1]
              if (dose %in% data@x){
                xd<-table(data@x==dose)
                if (dim(xd)>0){
                  n<-xd["TRUE"]
                  
                  #Load data events vector at this dose and find the number of DLTs
                  yd<-table(data@y,data@x==dose)
                  
                  #check if both DLTs and no DLTs are present
                  if (nrow(yd)==2){
                    #Number of DLTs
                    y1<-yd[2,"TRUE"]
                  }
                  #check if only DLTs or only no DLTs are present
                  if (nrow(yd)==1){
                    #Number of DLTs in each scenario:
                    #if only DLTs
                    if (dimnames(yd)[1]==1){y1=xd}
                    #if only no DLTs
                    if (dimnames(yd)[1]==0){y1=0}
                  }
                  #Calculate the probability of DLTs
                  # i.e. 1- pbeta(target, y+a, n-y+b) < 0.90
                  # where y=DLTs, n=patients for a dose and use Beta(a,b)
                  probTarget<-(1- pbeta(stopping@target, y1+stopping@a, n-y1+stopping@b))
                  
                  ## so can we stop?
                  doStop <- probTarget >= stopping@prob
                  
                  ## generate message
                  text <-
                    paste("Probability that the lowest dose",dose,"being toxic based on Beta(",
                          stopping@a,",",stopping@b,") is",probTarget,
                          "and thus",
                          ifelse(doStop, "above", "below"),
                          "the prespecified risk of",
                          stopping@prob)
                  
                  ## return both
                  return(structure(doStop,
                                   message=text))
                  
                }
              }else{
                
                probTarget<-NA
                
                ## so can we stop?
                doStop <- is.numeric(probTarget)
                
                ## generate message
                text <-
                  paste("Lowest dose",dose,"never tested before and thus probability the lowest dose being toxic based on Beta(",
                        stopping@a,",",stopping@b,") is",probTarget
                  )
                
                ## return both
                return(structure(doStop,
                                 message=text))
              }
              
              
            })
