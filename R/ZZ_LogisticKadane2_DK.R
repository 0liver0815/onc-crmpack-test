
## -----------------------------------------------------------------------------------------
## Add p0 , MTD in globalVariables in the program crmPack-package.r
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## New addition for [Model-class.R] 
## -----------------------------------------------------------------------------------------

## ============================================================

##' Reparametrized logistic model (version 2)
##'
##' This is the logistic model in the parametrization (version 2) of 
##' Kadane et al. (1980).
##'
##' Let \eqn{\rho_{0} = p(x_{0})} be the probability of a DLT of the placebo 
##' (no drug) dose \eqn{x_{0}}, and let \eqn{MTD} be the dose with target 
##' toxicity probability \eqn{\theta}, i.e. \eqn{p(MTD) = \theta}. Then it can 
##' easily be shown that the logistic regression model has 
##' intercept
##' \deqn{\frac{MTD logit(\rho_{0})}{MTD}}
##' and slope
##' \deqn{\frac{logit(\theta) - logit(\rho_{0})}{MTD}}
##' 
##' The prior for \eqn{MTD} is a \eqn{Gamma(shape,rate)} distribution.
##' The prior for \eqn{\rho_{0}} is a \eqn{Beta(\alpha,\beta)} distribution.
##' 
##' The minimum \eqn{d_{min}} and maximum \eqn{d_{max}} planned dose, are used 
##' to set the initial value of the \eqn{MTD} arbitrarily as the average of 
##' those two.
##' The initial value of \eqn{\rho_{0}} is set arbitrarily as \deqn{\frac{\theta}{10}}.
##'
##' The slots of this class, required for creating the model, are the target
##' toxicity, the Beta and Gamma distribution parameters, as well as the minimum
##' and maximum of the dose range. Note that these can be different from the 
##' minimum and maximum of the dose grid in the data later on.
##'
##' @slot theta the target toxicity probability \eqn{\theta}
##' @slot dmin the minimum of the dose range \eqn{d_{min}}
##' @slot dmax the maximum of the dose range \eqn{d_{max}}
##' @slot alpha the \eqn{\alpha} shape parameter of the \eqn{Beta(\alpha,\beta)} 
##' distribution
##' @slot beta the \eqn{\beta} shape parameter of the \eqn{Beta(\alpha,\beta)} 
##' distribution 
##' @slot shape the shape parameter of the \eqn{Gamma(shape,rate)} distribution 
##' @slot rate the rate parameter of the  \eqn{Gamma(shape,rate)} distribution.
##'
##' @export
##' @keywords classes

.LogisticKadane2 <- setClass(Class = "LogisticKadane2", 
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

validObject(.LogisticKadane2())

##' Initialization function for the "LogisticKadane2" class
##'
##' @param theta the target toxicity probability
##' @param dmin the minimum of the dose range 
##' @param dmax the maximum of the dose range 
##' @param alpha the \eqn{\alpha} shape parameter of the \eqn{Beta(\alpha,\beta)} 
##' distribution
##' @param beta the \eqn{\beta} shape parameter of the \eqn{Beta(\alpha,\beta)} 
##' distribution 
##' @param shape the shape parameter of the \eqn{Gamma(shape,rate)} distribution 
##' @param rate the rate parameter of the  \eqn{Gamma(shape,rate)} distribution.
##' @return the \code{\linkS4class{LogisticKadane2}}
##'
##' @export
##' @keywords methods
##' 
LogisticKadane2 <- function(theta, dmin, dmax, alpha, beta, shape, rate)
{
  .LogisticKadane2(theta = theta, 
                   dmin = dmin, 
                   dmax = dmax, 
                   alpha = alpha, 
                   beta = beta, 
                   shape = shape, 
                   rate = rate,
                   datamodel = function() {
                     for (i in 1:nObs) {
                       y[i] ~ dbern(p[i])
                       logit(p[i]) <- (1/(MTD)) * (MTD * logit(p0) + 
                                                     (logit(theta) - 
                                                        logit(p0)) * x[i])
                       ## the logit(p[i]) definition above is the same as the form
                       #  logit(p[i]) <- logit(p0) + ((logit(theta)-logit(p0))/MTD)*x[i]
                       #  used in the protocol appendix
                     }
                   }, 
                   priormodel = function() {
                     p0 ~ dbeta(alpha, beta)
                     MTD ~ dgamma(shape, rate)
                     
                     ## the dummy assignment below is arbitrary and is just added to eliminate warning messages input variables not being used anywhere in the function
                     #  Warning messages:
                     #   1: In rjags::jags.model(file = modelFileName, data = if (fromPrior) modelspecs else c(requiredData,  :
                     #      Unused variable "dmin" in data
                     #   2: In rjags::jags.model(file = modelFileName, data = if (fromPrior) modelspecs else c(requiredData,  :
                     #      Unused variable "dmax" in data
                     #   3. In rjags::jags.model(file = modelFileName, data = if (fromPrior) modelspecs else c(requiredData,  :
                     #      Unused variable "theta" in data
                     
                     lowestdose <- dmin
                     highestdose <- dmax
                     DLTtarget <- theta
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
                     ret <- MTD * (logit(prob) - logit(p0))
                     ret <- ret/(logit(theta) - logit(p0))
                     return(ret)
                     ## the definition above is the same as
                     #  (MTD * (logit(prob) - logit(p0))) / (logit(theta) - logit(p0))
                   }, 
                   prob = function(dose, p0, MTD) {
                     ret <- (MTD * logit(p0) + 
                               (logit(theta) - logit(p0)) * dose)
                     ret <- plogis(ret/MTD)
                     return(ret)
                     ## the definition above is the same as
                     #  1/(1 + exp(-logit(p0) - ((logit(theta)-logit(p0))/MTD) * dose))
                     #  plogis is the R logistic distribution function
                     
                   }, 
                   init = function() {
                     list(p0 = theta/10, MTD = (dmax + dmin)/2)
                     # init = function() {
                     #  list(p0 = theta/10, MTD = (max(data@doseGrid[]) - data@doseGrid[1])/2)
                   }, 
                   sample = c("p0", "MTD"))
}