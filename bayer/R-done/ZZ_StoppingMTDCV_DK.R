## -----------------------------------------------------------------------------------------
## Add mad in importFrom stats in the program crmPack-package.r
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## New addition for [Rules-class.R] 
## -----------------------------------------------------------------------------------------

## --------------------------------------------------
## Stopping based on precision of MTD calculated as CV(MTD)
## --------------------------------------------------

##' Stop based on precision of MTD calculated as CV(MTD)
##'
##' @slot target toxicity target of MTD
##' @slot thresh threshold for MTD to be considered accurate enough and stop 
##' the trial
##' 
##' @keywords classes
##' @export

.StoppingMTDCV_dk <-
  setClass(Class="StoppingMTDCV_dk",
           representation(target="numeric",
                          thresh="numeric"),
           prototype(target=0.33,
                     thresh=0.4),
           contains="Stopping",
           validity=
             function(object){
               o <- Validate()
               
               o$check(is.probability(object@target,
                                      bounds=FALSE),
                       "target must be probability > 0 and < 1")
               o$check(is.probability(object@thresh,
                                      bounds=FALSE),
                       "thresh must be probability > 0 and < 1")
               
               o$result()
             })
validObject(.StoppingMTDCV_dk())

##' Initialization function for "StoppingMTDCV_dk"
##'
##' @param target see \code{\linkS4class{StoppingMTDCV_dk}}
##' @param thresh see \code{\linkS4class{StoppingMTDCV_dk}}
##' @return the \code{\linkS4class{StoppingMTDCV_dk}} object
##'
##' @export
##' @keywords methods

StoppingMTDCV_dk <- function(target,
                          thresh)
{
  .StoppingMTDCV_dk(target=target,
                 thresh=thresh)
}

## -----------------------------------------------------------------------------------------
## New addition for [Rules-methods.R] 
## -----------------------------------------------------------------------------------------

## --------------------------------------------------
## Stopping based on precision of MTD calculated as CV(MTD)
## --------------------------------------------------

##' @describeIn stopTrial Stop based on precision of MTD calculated as CV(MTD)

setMethod("stopTrial",
          signature=
            signature(stopping="StoppingMTDCV_dk",
                      dose="numeric",
                      samples="Samples",
                      model="Model",
                      data="ANY"),
          def=
            function(stopping, dose, samples, model, data, ...){
              ## First, generate the MTD samples.
              
              ## add prior data and samples to the
              ## function environment so that they
              ## can be used.
              MTDSamples <- dose(prob=stopping@target,
                                 model,
                                 samples)
              
              ## CV of MTD derived based on MTD samples
              MTDCV <- (mad(MTDSamples)/median(MTDSamples))
              
              ## so can we stop?
              doStop <- ((MTDCV <= stopping@thresh) & (MTDCV >= 0))
              
              ## generate message
              text <-
                paste("CV of MTD is",
                      round(MTDCV * 100),"% and thus",
                      ifelse(doStop, "below", "above"),
                      "the required precision threshold of",
                      round(stopping@thresh * 100),
                      "%")
              
              ## return both
              return(structure(doStop,
                               message=text))
            })
