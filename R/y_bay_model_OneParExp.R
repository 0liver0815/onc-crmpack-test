#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# One parameter model: model and MTD selection rule 
#
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------------------------------------------------------------------------------
# Implementation of the one parameter power model
#-------------------------------------------------------------------------------

##' Standard one parameter model
##'
##' @slot skeletonFun the skeleton
##' @slot skeletonProbs the skeleton probs
##' @slot lambda the parameter lambda
##'
##' @export
##' @keywords classes
.OneParExp <- setClass(Class = "OneParExp",
                       contains = "Model",
                       representation(skeletonFun = "function",
                                      skeletonProbs = "numeric",
                                      lambda = "numeric"))

##' Initialization function for the "OneParExp" class
##'
##' @param skeletonProbs the skeleton probs
##' @param doseGrid the dose grid
##' @param lambda the lambda
##' @return the \code{\linkS4class{OneParExp}} object
##'
##' @importFrom stats approxfun
##' @export
##' @keywords methods
OneParExp <- function(skeletonProbs, doseGrid, lambda) {
  skeletonFun <- approxfun(x = doseGrid, y = skeletonProbs, rule = 2)
  invSkeletonFun <- approxfun(x = skeletonProbs, y = doseGrid, rule = 1)
  .OneParExp(skeletonFun = skeletonFun, skeletonProbs = skeletonProbs, lambda = lambda,
             datamodel = function() {
               for (i in 1:nObs) {
                 y[i] ~ dbern(p[i])
                 p[i] <- skeletonProbs[xLevel[i]]^theta
               }
             },
             datanames = c("nObs", "y", "xLevel"),
             prob = function(dose, theta) {
               skeletonFun(dose)^theta
             },
             dose = function(prob, theta) {
               invSkeletonFun(prob^(1 / theta))
             },
             priormodel = function() {
               theta ~ dexp(lambda)
             },
             modelspecs = function() {
               list(skeletonProbs = skeletonProbs, lambda = lambda)
             },
             init = function() {
               list(theta = 1)
             }, sample = "theta")
}
