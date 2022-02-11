#' Stopping rule based on the posterior distribution. The trial is stopped,
#' when the MTD can be estimated with sufficient precision.
#' The criteria is based on the robust CV calculated from the posterior
#' distribution. The robust coefficient of variation used is defined
#' as MAD(MTD)/median(MTD)
#'
#' @slot target the target toxicity
#' @slot CVthreshold the minimum precision
#' @slot label label of stopping rule
#' 
#' @export
#' @keywords classes
.StoppingCV <-
  setClass(Class="StoppingCV",
           representation(target="numeric",
                          CVthreshold="numeric",
                          label="character"),
           prototype(target=0.33,
                     CVthreshold=40,
                     label="precision (CV)"),
           contains="Stopping")


#' Initialization function for "StoppingCV"
#'
#' @param target see \code{\linkS4class{StoppingCV}}
#' @param CVthreshold see \code{\linkS4class{StoppingCV}}
#' @param label see \code{\linkS4class{StoppingCV}}
#'
#' @export
#' @keywords methods
StoppingCV <- function(target,
                       CVthreshold,
                       label)
{
  .StoppingCV(target=target,
              CVthreshold=CVthreshold,
              label=label)
}

#' @describeIn stopTrial Stopping rule based precision of the MTD estimation
#' The trial is stopped, when the MTD can be estimated with sufficient precision.
#' The criteria is based on the robust CV calculated from the posterior distribution.
#' The robust coefficient of variation is defined as MAD(MTD)/median(MTD)


#' @importFrom stats mad
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingCV",
                      dose="numeric",
                      samples="Samples"),
          def=
            function(stopping, dose, samples, model, data, ...){
              
              # add prior data and samples to the
              # function environment so that they
              # can be used.
              mtdSamples <- dose(prob=stopping@target,
                                 model,
                                 samples)
              # calculate robust coefficient of variation
              #mtd <- (log(stopping@target/(1-stopping@target)) - samples@data$int) / samples@data$slope
              #cv <- mad(mtd)/median(mtd)
              cv <- mad(mtdSamples)/median(mtdSamples)
              
              # so can we stop?
              doStop <- cv*100 <= stopping@CVthreshold
              
              # generate message
              text <-
                paste(doStop, ' : ', stopping@label,
                      " : Robust coefficient of variation (MAD/Median) is ",
                      round(cv * 100),
                      "% and thus ",
                      ifelse(doStop, "below", "above"),
                      " the required ",
                      round(stopping@CVthreshold),
                      "%",
                      sep="")
              
              # return both
              return(structure(doStop,
                               message=text))
            })
