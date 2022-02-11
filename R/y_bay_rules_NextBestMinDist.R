#-------------------------------------------------------------------------------
# Implementation of the selection rule: minimal distance to target
#-------------------------------------------------------------------------------

#' The class with the input for finding the next best MTD estimate with minimal distance
#'
#' @slot target the target toxicity probability
#' 
#' @export
#' @keywords classes
.NextBestMinDist <- setClass(Class = "NextBestMinDist",
                             contains = "NextBest",
                             representation(target = "numeric"))


#' Initialization function for class "NextBestMinDist"
#'
#' @param target see \code{\linkS4class{NextBestMTD}}
#'
#' @export
#' @keywords methods
NextBestMinDist <- function(target) {
  .NextBestMinDist(target = target)
}

#' @describeIn nextBest Find the next best dose based on the minimal distance to
#' the target
setMethod("nextBest",
          signature = signature(nextBest = "NextBestMinDist",
                                doselimit = "numeric",
                                samples = "Samples",
                                model = "Model",
                                data = "Data"),
          def = function(nextBest, doselimit, samples, model, data, ...) {
            dosesOK <- if (length(doselimit)) {
              which(data@doseGrid <= doselimit)
            } else {
              seq_along(data@doseGrid)
            }
            modelfit <- fit(samples, model, data)
            probDLT <- modelfit$middle[dosesOK]
            doses <- modelfit$dose[dosesOK]
            bestIndex <- which.min(abs(probDLT - nextBest@target))
            bestDose <- doses[bestIndex]
            return(list(value = bestDose))
          })

