#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# Two parameter model Bayer implementation: model and MTD selection rule 
#
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------------------------------------------------------------------------------
# Implementation of the two parameter logistic model without dose scaling
#-------------------------------------------------------------------------------

##' Standard two parameter model
##'
##' @slot mean1 the mean of the intercept
##' @slot mean2 the mean of the slope
##' @slot var1 the variance of the intercept
##' @slot var2 the variance of the slope
##'
##' @export
##' @keywords classes
.TwoParBay <- setClass(Class = "TwoParBay",
                       contains = "Model",
                       representation(mean1 = "numeric",
                                      mean2 = "numeric",
                                      var1 = "numeric",
                                      var2 = "numeric"))


globalVariables(c("int",
                  "slope")) 

##' Initialization function for the "TwoParBay" class
##'
##' @param mean1 the mean of the intercept
##' @param mean2 the mean of the slope
##' @param var1 the variance of the intercept
##' @param var2 the variance of the slope
##' @return the \code{\linkS4class{TwoParBay}} object
##'
##' @export
##' @keywords methods
TwoParBay <- function (mean1,
                       mean2,
                       var1,
                       var2) 
{
  .TwoParBay( mean1 = mean1,
              mean2 = mean2,
              var1  = var1,
              var2  = var2,
              datamodel = 
                function() {
                  for (i in 1:nObs) 
                  {
                    y[i] ~ dbern(mean[i])
                    logit(mean[i]) <- int + slope * x[i]
                  }
                },
              priormodel =
                function() {
                  int ~ dnorm (mean1, var1)
                  slope ~ dnorm (mean2, var2)%_%I(0,)
                },
              datanames = c("nObs", "y", "x"), 
              modelspecs =
                function() {
                  list(mean1 = mean1,
                       mean2 = mean2,
                       var1  = var1,
                       var2  = var2)
                },
              dose =
                function(prob, int, slope) {
                  (logit(prob)-int)/slope
                },
              prob =
                function(dose, int, slope) {
                  1/(1+exp(-int-slope*dose))
                },
              init = 
                function() {
                  list(int = mean1, slope = mean2)
                },
              sample = c("int", "slope"))
}
