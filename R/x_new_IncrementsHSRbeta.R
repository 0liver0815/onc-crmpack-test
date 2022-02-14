#
#-->rules-validity.R
#
#' @describeIn validate_stopping validates that the [`IncrementsHSRBeta`]
#'  object contains valid probability target, threshold and shape parameters.
validate_increments_hrs_beta <- function(object) {
  o <- Validate()
  o$check(
    is.probability(object@target, bounds = FALSE),
    "target must be probability > 0 and < 1"
  )
  o$check(
    is.probability(object@prob, bounds = FALSE),
    "prob must be probability > 0 and < 1"
  )
  o$check(
    is.numeric(object@a) & object@a > 0,
    "Beta distribution shape parameter a must me a real number > 0"
  )
  o$check(
    is.numeric(object@b) & object@b > 0,
    "Beta distribution shape parameter b must me a real number > 0"
  )

  o$result()
}

# IncrementsHSRBeta-class ----

#' `IncrementsHSRBeta`
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' [`IncrementsHSRBeta`] is a class for increments based on the Bin-Beta model.
#' Increment control is based on the number of observed DLTs and number of
#' subjects at each dose level. The probability of toxicity is calculated
#' using a Bin-Beta model with prior (a,b). If the probability exceeds
#' the threshold for a given dose, that dose and all doses above are excluded
#' from further escalation.
#' This is a hard safety rule that limits further escalation based on the
#' observed data per dose level.
#'
#' @slot target (`number`)\cr the target toxicity
#' @slot prob (`number`)\cr the threshold probability for a dose being toxic
#' @slot a (`number`)\cr shape parameter a>0 of probability
#'  distribution Beta (a,b)
#' @slot b (`number`)\cr shape parameter b>0 of probability
#'  distribution Beta (a,b)
#'
#' @aliases IncrementsHSRBeta
#' @export
#'
.IncrementsHSRBeta <- setClass(
  Class = "IncrementsHSRBeta",
  contains = "Increments",
  representation(
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
  validity = validate_increments_hrs_beta
)

# IncrementsHSRBeta-constructor ----

#' @rdname IncrementsHSRBeta-class
#'
#' @param target (`number`)\cr the target toxicity
#' @param prob (`number`)\cr the threshold probability for a dose being toxic
#' @param a (`number`)\cr shape parameter a>0 of probability
#'  distribution Beta (a,b)
#' @param b (`number`)\cr shape parameter b>0 of probability
#'  distribution Beta (a,b)
#'
#' @example bayer/examples/Rules-class-IncrementsHSRBeta.R
#' @export
#'
IncrementsHSRBeta <- function(target = 0.3,
                              prob = 0.95,
                              a = 1,
                              b = 1) {
  .IncrementsHSRBeta(
    target = target,
    prob = prob,
    a = a,
    b = b
  )
}


# maxDose-IncrementsHSRBeta ----

#' @rdname maxDose
#'
#' @description Determine the maximum possible dose for escalation.
#'
#' @aliases maxDose-IncrementsHSRBeta
#' @example bayer/examples/Rules-method-maxDose-IncrementsHSRBeta.R
#' @export
setMethod(
  "maxDose",
  signature = signature(
    increments = "IncrementsHSRBeta",
    data = "Data"
  ),
  definition = function(increments, data, ...) {
    # summary of observed data per dose level
    y <- factor(data@y, levels = c("0", "1"))
    dlt_tab <- table(y, data@x)

    # ignore placebo if applied
    if (data@placebo==TRUE & min(data@x) == data@doseGrid[1]){
      dlt_tab <- dlt_tab[,-1]
    }

    # extract dose names as these get lost if only one dose available
    doses <- as.numeric(colnames(dlt_tab))

    # toxicity probability per dose level
    tox_prob <- 1 - pbeta(
      increments@target, dlt_tab[2, ] + increments@a,
      apply(dlt_tab, 2, sum) - dlt_tab[2, ] + increments@b
    )

    # return the min toxic dose level or maximum dose level if no dose is toxic
    # while ignoring placebo
    if (sum(tox_prob > increments@prob) > 0) {
      dose_tox <- min(doses[which(tox_prob > increments@prob)])
    } else {
      # add small value to max dose, so that the max dose is always smaller
      dose_tox <- max(data@doseGrid) + 0.01
    }

    # Determine the next maximum possible dose.
    # In case that the first active dose is above probability threshold,
    # the first active dose is reported as maximum. I.e. in case that placebo is used,
    # the second dose is reported. Please note that this rule should be used together
    # with the hard safety stopping rule to avoid inconsistent results.
    max(data@doseGrid[data@doseGrid < dose_tox], data@doseGrid[data@placebo+1])
  }
)
