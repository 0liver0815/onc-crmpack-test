#' Stop based on precision of MTD calculated as CV(MTD)
#' 
#' @description `r lifecycle::badge("experimental")`
#'
#' @slot target (`numeric`)\cr toxicity target of MTD
#' @slot threshCV (`numeric`)\cr threshold for CV to be considered accurate enough
#'  to stop the trial
#' 
#' @example examples_new/Rules-class-StoppingMTDCV.R
#' @export
.StoppingMTDCV <-
  setClass(
    Class = "StoppingMTDCV",
    representation(
      target = "numeric",
      threshCV = "numeric"
    ),
    prototype(
      target = 0.33,
      threshCV = 40
    ),
    contains = "Stopping",
    validity =
      function(object) {
        o <- Validate()
        o$check(
          is.probability(object@target, bounds = FALSE),
          "target must be probability > 0 and < 1"
        )
        o$check(
          is.probability(object@threshCV/100, bounds = FALSE),
          "threshCV must be percentage > 0"
        )
        o$result()
      }
  )
validObject(.StoppingMTDCV())

#' Initialization function for `StoppingMTDCV`
#'
#' @param target see above.
#' @param threshCV see above.
#' @return the `StoppingMTDCV` object.
#'
#' @rdname StoppingMTDCV-class
#' @export
StoppingMTDCV <- function(target,
                          threshCV) {
  .StoppingMTDCV(
    target = target,
    threshCV = threshCV
  )
}

## --------------------------------------------------
## Stopping based on precision of MTD calculated as CV(MTD)
## --------------------------------------------------

#' @describeIn stopTrial Stopping rule based precision of the MTD estimation.
#' The trial is stopped, when the MTD can be estimated with sufficient precision.
#' The criteria is based on the robust CV calculated from the posterior distribution.
#' The robust coefficient of variation is defined as MAD(MTD)/median(MTD).
#' @example examples_new/Rules-method-stopTrial-StoppingMTDCV.R
setMethod("stopTrial", 
          signature = signature(
            stopping = "StoppingMTDCV",
            dose = "numeric",
            samples = "Samples",
            model = "Model",
            data = "ANY"
          ),
          def =
            function(stopping, dose, samples, model, data, ...) {
              # First, generate the MTD samples.
              # add prior data and samples to the
              # function environment so that they
              # can be used.
              mtd_samples <- dose(
                prob = stopping@target,
                model,
                samples
              )
              
              # CV of MTD expressed as percentage, derived based on MTD posterior samples
              mtd_cv <- (mad(mtd_samples) / median(mtd_samples)) * 100
              
              # so can we stop?
              do_stop <- ((mtd_cv <= stopping@threshCV) & (mtd_cv >= 0))
              
              # generate message
              msg <- paste(
                "CV of MTD is ",
                round(mtd_cv), "% and thus ",
                ifelse(do_stop, "below", "above"),
                " the required precision threshold of ",
                round(stopping@threshCV),
                "%",
                sep=''
              )
              
              # return both
              structure(do_stop, message = msg)
            }
)

