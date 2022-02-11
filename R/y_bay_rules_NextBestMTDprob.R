#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Next best dose
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' The class with the input for finding the next best MTD estimate with highest
#' probability to be the MTD.
#' The dose which has the highest probability to be the MTD among the specified
#' doses is used as next best dose. The posterior is calculated for each
#' iteration and the maximum dose that is not toxic (below target) is determined.
#' The allocation criteria is calculated as the frequency of a dose being
#' selected as MTD. The dose with the highest allocation value is selected as
#' next best dose.
#' 
#' @slot target the target toxicity probability
#' @slot method the method used for calculating the allocation criteria
#' 1=aggregate no dose safe with first dose
#' 2=aggregate no dose safe with last dose
#' 3=do not aggregate, but report the first dose as best next dose, in case 'no
#'   dose is safe' has the highest allocation value. This follows the idea used for 
#'   the smallest distance to the target dose approach, where always a dose is
#'   recommended, even though the lowest dose estimated is above target
#' 
#' @export
#' @keywords classes
.NextBestMTDprob <- setClass(Class = "NextBestMTDprob",
                             contains = "NextBest",
                             representation(target = "numeric",
                                            method="numeric"),
                             prototype(target = 0.3,
                                       method = 1),
)


#' Initialization function for class "NextBestMTDprob"
#'
#' @param target see \code{\linkS4class{NextBestMTDprob}}
#' @param method see \code{\linkS4class{NextBestMTDprob}}
#'
#' @export
#' @keywords methods
NextBestMTDprob <- function(target,
                            method=1) {
  .NextBestMTDprob(target = target,
                   method = method)
}

#' @describeIn nextBest Find the next best dose based on MTD estimate with highest
#' probability to be the MTD
setMethod("nextBest",
          signature = signature(nextBest = "NextBestMTDprob",
                                doselimit = "numeric", samples = "Samples", model = "Model",
                                data = "Data"),
          def = function(nextBest, doselimit, samples, model, data, ...) {
            
            # be sure which doses are ok with respect to maximum
            # possible dose - if one was specified
            dosesOK <-
              if(length(doselimit))
                which(data@doseGrid <= doselimit)
            else
              seq_along(data@doseGrid)
            
            # determine maximum allowed index
            # maxdoseOK <- if (length(doselimit)) {
            #   max(which(data@doseGrid <= doselimit),0)
            # } else {data@nGrid}
            
            # Number of iterations from mcmc sampling
            sampleLength <- (samples@options@iterations-samples@options@burnin) / samples@options@step
            
            
            # create an empty matrix of the dimension nGrid x sampleLength
            probdose <- matrix(NA,
                               ncol=sampleLength,
                               nrow=data@nGrid)
            
            # extract the probability function from the model 
            probFun <- slot(model, "prob")
            
            # which arguments, besides the dose, does it need?
            argNames <- setdiff(names(formals(probFun)), "dose")
            
            # calculate the toxicity for any sample and dose:
            # now call the function with any dose in the doseGrid and with
            # the arguments taken from the samples and store the results
            # in the columns of the matrix
            for (i in seq_along(data@doseGrid)){
              probdose[i,] <- do.call(probFun,
                                      c(list(dose=data@doseGrid[i]),
                                        samples@data[argNames]))
            }
            
            # determine which dose level is just below the target
            # and count the rel frequency of each dose level
            probMTDdist <- table(factor(apply(probdose < nextBest@target, 2, sum),
                                        levels = 0:data@nGrid)) / sampleLength
            
            # aggregate the relative frequency according to the selected method
            # If method=1, aggregate no dose below target with dose 1 below target
            # if method=2 aggregate no dose below target and last dose below target
            # if method=3 do not aggregate, but take same decision as if first dose
            if (nextBest@method == 1){
              probMTDdist[2] <- sum(probMTDdist[1:2])
            }else if (nextBest@method == 2){
              probMTDdist[data@nGrid+1] <- sum(probMTDdist[c(1, data@nGrid+1)])
            }else if (nextBest@method == 3){
              probMTDdist[2] <- max(probMTDdist[1:2])
            }
            # remove no dose is below target
            probMTDdist <- probMTDdist[2:length(probMTDdist)]
            
            # determine the dose with the highest frequency and is allowed
            # if there is no doseOK the following results can be triggered:
            # max(dosesOK,-Inf) leads to NA as nextbest, but than no stop
            # rule is hit
            # max(dosesOK,1) leads to selection of dose 1 and stopping rules
            # will be evaluated
            bestIndex <- as.integer(names(which.max(probMTDdist)))
            bestIndex <- min(bestIndex, max(dosesOK,1))
            
            # determine the next best dose
            bestDose <- data@doseGrid[bestIndex]
            
            # replace index by actual dose values
            names(probMTDdist) <- data@doseGrid 
            
            return(list(value = bestDose,
                        allocation=probMTDdist,
                        index=bestIndex))
          })

