#
#-->rules-validity.R
#
#' @describeIn validate_stopping validates that the [`StoppingLowestDoseHSRBeta`]
#'  object contains valid probability target, threshold and shape parameters.
validate_stopping_lowest_dose_hsr_beta <- validate_increments_hrs_beta
# validate_stopping_lowest_dose_hsr_beta <- function(object) {
#   o <- Validate()
#   o$check(
#     is.probability(object@target, bounds = FALSE),
#     "target must be probability > 0 and < 1"
#   )
#   o$check(
#     is.probability(object@prob, bounds = FALSE),
#     "prob must be probability > 0 and < 1"
#   )
#   o$check(
#     is.numeric(object@a) & object@a > 0,
#     "Beta distribution shape parameter a must me a real number > 0"
#   )
#   o$check(
#     is.numeric(object@b) & object@b > 0,
#     "Beta distribution shape parameter b must me a real number > 0"
#   )
#   o$result()
# }


# StoppingLowestDoseHSRBeta-class ----

#' `StoppingLowestDoseHSRBeta`
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' [`StoppingLowestDoseHSRBeta`] is a class for stopping based on a hard safety
#' rule using the Beta posterior distribution with Beta(a,b) prior and a
#' Bin-Beta model based on the observed data at the lowest dose level.
#' The rule is triggered when the first dose is considered to be toxic
#' (i.e. above threshold probability) based on the observed data at the
#' lowest dose level and a Beta(a,b) prior distribution.
#' The default prior is Beta(1,1).
#' In case that placebo is used, the rule is evaluated at the second dose of the
#' dosegrid, i.e. at the lowest non-placebo dose.
#'
#' @slot target (`number`)\cr the target toxicity
#' @slot prob (`number`)\cr the threshold probability for the lowest
#'  dose being toxic
#' @slot a (`number`)\cr shape parameter a>0 of probability
#'  distribution Beta (a,b)
#' @slot b (`number`)\cr shape parameter b>0 of probability
#'  distribution Beta (a,b)
#'
#' @aliases StoppingLowestDoseHSRBeta
#' @export
#'
.StoppingLowestDoseHSRBeta <- setClass(
  Class = "StoppingLowestDoseHSRBeta",
  contains = "Stopping",
  representation = representation(
    target = "numeric",
    prob = "numeric",
    a = "numeric",
    b = "numeric"
  ),
  prototype(
    target = 0.3,
    prob = 0.95,
    a = 1,
    b = 1
  ),
  validity = validate_stopping_lowest_dose_hsr_beta
)


# StoppingLowestDoseHSRBeta-constructor ----

#' @rdname StoppingLowestDoseHSRBeta-class
#'
#' @param target (`number`)\cr the target toxicity probability (e.g. `0.3`)
#'   defining the MTD.
#' @param prob (`number`)\cr the threshold probability for the lowest
#'  dose being toxic (e.g. `0.95`).
#' @param a (`number`)\cr shape parameter a>0 of probability
#'  distribution Beta (a,b), default is Beta(1,1)
#' @param b (`number`)\cr shape parameter b>0 of probability
#'  distribution Beta (a,b), default is Beta(1,1)
#'
#' @export
#' @example bayer/examples/Rules-class-StoppingLowestDoseHSRBeta.R
#'
StoppingLowestDoseHSRBeta <- function(target = 0.3,
                                      prob = 0.95,
                                      a = 1,
                                      b = 1) {
  .StoppingLowestDoseHSRBeta(
    target = target,
    prob = prob,
    a = a,
    b = b
  )
}

# stopTrial-StoppingLowestDoseHSRBeta ----

#' @rdname stopTrial
#'
#' @description Stopping based based on the lowest dose meeting the hard
#' safety criteria using Bin-Beta model based DLT probability.
#'
#' @aliases stopTrial-StoppingLowestDoseHSRBeta
#' @example bayer/examples/Rules-method-stopTrial-StoppingLowestDoseHSRBeta.R
#' @export
setMethod(
  "stopTrial",
  signature = signature(
    stopping = "StoppingLowestDoseHSRBeta",
    dose = "numeric",
    samples = "Samples"
  ),
  definition = function(stopping, dose, samples, model, data, ...) {
    # determine if the first doses is toxic
    # First active dose Tested?
    if (sum(data@x == data@doseGrid[data@placebo+1]) > 0) {
      lowest_dose_tested <- TRUE
      # summary data at first dose
      y <- factor(data@y, levels = c("0", "1"))
      dlt_tab <- table(y, data@x)[, data@placebo+1]
      tox_prob_first_dose <-
          1 - pbeta(
            stopping@target, dlt_tab[2] + stopping@a,
            sum(dlt_tab) - dlt_tab[2] + stopping@b
          )
    } else {
      lowest_dose_tested <- FALSE
      tox_prob_first_dose <- 0
    }

    # so can we stop?
    do_stop <- tox_prob_first_dose > stopping@prob

    # generate message
    msg <- if (lowest_dose_tested == FALSE) {
      paste("Lowest active dose not tested, stopping rule not applied.")
    } else {
      paste(
        "Probability that the lowest active dose of ",
        data@doseGrid[data@placebo+1],
        " being toxic based on posterior Beta distribution using a Beta(",
        stopping@a, ",", stopping@b, ") prior is ",
        round(tox_prob_first_dose * 100),
        "% and thus ",
        ifelse(do_stop, "above", "below"),
        " the required ",
        round(stopping@prob * 100),
        "% threshold.",
        sep = ""
      )
    }

    # return both
    structure(do_stop, message = msg)
  }
)
