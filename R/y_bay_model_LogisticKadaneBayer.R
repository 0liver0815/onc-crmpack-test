#' Standard Kadane model
#'
#' @slot theta theta
#' @slot dmin dose min
#' @slot dmax dose max
#' @slot alpha intercept
#' @slot beta slope
#' @slot shape shape
#' @slot rate rate
#'
#' @export
#' @keywords classes
#' 
.LogisticKadaneBayer <- setClass(Class = "LogisticKadaneBayer",
                                 contains = "Model",
                                 representation(theta = "numeric",
                                                dmin = "numeric",
                                                dmax = "numeric",
                                                alpha = "numeric",
                                                beta = "numeric",
                                                shape = "numeric",
                                                rate = "numeric"),
                                 prototype(theta=0.3,
                                           dmin=0.1,
                                           dmax=1,
                                           alpha=1,
                                           beta=0.5,
                                           shape=1.2,
                                           rate=2.5),
                                 validity=
                                   function(object){
                                     o <- Validate()
                                     
                                     o$check(is.probability(object@theta,
                                                            bounds=FALSE),
                                             "theta must be a probability > 0 and < 1")
                                     o$check(object@dmin < object@dmax,
                                             "dmin must be smaller than dmax")
                                     o$check(is.scalar(object@dmin),
                                             "dmin must be scalar")
                                     o$check(is.scalar(object@dmax),
                                             "dmax must be scalar")
                                     o$check(object@alpha > 0,
                                             "alpha must be greater than 0")
                                     o$check(object@beta > 0,
                                             "beta must be greater than 0")
                                     o$check(object@shape > 0,
                                             "shape must be greater than 0")
                                     o$check(object@rate > 0,
                                             "rate must be greater than 0")
                                     o$result()
                                   })


globalVariables(c("p0", "MTD")) 

#' Initialization function for the "LogisticKadaneBayer" class
#'
#' @param theta theta
#' @param dmin dose min
#' @param dmax dose max
#' @param alpha intercept
#' @param beta slope
#' @param shape shape
#' @param rate rate
#' @return the \code{\linkS4class{LogisticKadaneBayer}} object
#'
#' @export
#' @keywords methods
LogisticKadaneBayer <- function(theta, dmin, dmax, alpha, beta, shape, rate)
{
  .LogisticKadaneBayer(theta = theta, 
                       dmin = dmin, 
                       dmax = dmax, 
                       alpha = alpha, 
                       beta = beta, 
                       shape = shape, 
                       rate = rate,
                       datamodel = function() {
                         for (i in 1:nObs) {
                           y[i] ~ dbern(p[i])
                           logit(p[i]) <- logit(p0) + ((logit(theta)-logit(p0))/MTD)*x[i] 
                         }
                       }, 
                       priormodel = function() {
                         p0 ~ dbeta(alpha, beta)
                         MTD ~ dgamma(shape, rate)
                         ## dummy to use dmin & dmax
                         ## It is contained in the modelspecs list below,
                         ## so it must occur here
                         bla <- dmin + dmax + theta
                       }, 
                       datanames = c("nObs", "y", "x"), 
                       modelspecs = function() {
                         list(theta = theta, 
                              dmin = dmin, 
                              dmax = dmax, 
                              alpha = alpha, 
                              beta = beta, 
                              shape = shape, 
                              rate = rate)
                       }, 
                       dose = function(prob, p0, MTD) {
                         (MTD * (logit(prob) - logit(p0))) / (logit(theta) - logit(p0))
                       }, 
                       prob = function(dose, p0, MTD) {
                         1/(1 + exp(-logit(p0) - ((logit(theta)-logit(p0))/MTD) * dose))
                       }, 
                       init = function() {
                         list(p0 = theta/10, MTD = (dmax + dmin)/2)
                       }, 
                       sample = c("p0", "MTD"))
}

