# The new code below is necessary as no MTD estimate is returned by the simulate
# function. The below code introduce a new class BayDesign, which is simulated with
# a different simulation function that return the MTD estimate as last row
# in the toxicity estimated from the fit function

# --------------------------------------------------
# Classes for model-based designs Bayer Version
# --------------------------------------------------

#' Class for the Bayer CRM design
#'
#' In addition to the slots in the more simple \code{\linkS4class{RuleDesign}},
#' objects of this class contain:
#'
#' @slot model the model to be used, an object of class
#' \code{\linkS4class{Model}}
#' @slot stopping stopping rule(s) for the trial, an object of class
#' \code{\linkS4class{Stopping}}
#' @slot increments how to control increments between dose levels,
#' an object of class \code{\linkS4class{Increments}}
#' @slot PLcohortSize rules for the cohort sizes for placebo, if any planned
#' an object of class \code{\linkS4class{CohortSize}} (defaults to constant
#' 0 placebo patients)
#'
#' @example examples/design-class-Design.R
#' @export
#' @keywords classes
.BayDesign <-
  setClass(Class="BayDesign",
           representation(model="Model",
                          stopping="Stopping",
                          increments="Increments",
                          PLcohortSize="CohortSize"),
           prototype(model=.LogisticNormal(),
                     nextBest=.NextBestNCRM(),
                     stopping=.StoppingMinPatients(),
                     increments=.IncrementsRelative(),
                     PLcohortSize=CohortSizeConst(0L)),
           contains=list("RuleDesign"))
validObject(.Design())


#' Initialization function for "Design"
#'
#' @param model see \code{\linkS4class{Design}}
#' @param stopping see \code{\linkS4class{Design}}
#' @param increments see \code{\linkS4class{Design}}
#' @param PLcohortSize see \code{\linkS4class{Design}}
#' @param \dots additional arguments for \code{\link{RuleDesign}}
#' @return the \code{\linkS4class{BayDesign}} object
#'
#' @export
#' @keywords methods
BayDesign <- function(model,
                      stopping,
                      increments,
                      PLcohortSize=CohortSizeConst(0L),
                      ...)
{
  start <- RuleDesign(...)
  .BayDesign(start,
             model=model,
             stopping=stopping,
             increments=increments,
             PLcohortSize=PLcohortSize)
}


