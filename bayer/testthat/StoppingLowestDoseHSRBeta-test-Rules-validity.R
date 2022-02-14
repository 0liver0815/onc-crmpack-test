# validate_stopping_lowest_dose_hsr_beta ----

test_that("validate_stopping_lowest_dose_hsr_beta passes for valid object", {
  object <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.95)
  expect_true(validate_stopping_lowest_dose_hsr_beta(object))
})

test_that("validate_stopping_lowest_dose_hsr_beta passes for valid object", {
  object <- StoppingLowestDoseHSRBeta(target = 0.2, prob = 0.9, a = 7, b = 3)
  expect_true(validate_stopping_lowest_dose_hsr_beta(object))
})

test_that("StoppingLowestDoseHSRBeta returns expected messages for non-valid object", {
  object <- StoppingLowestDoseHSRBeta()
  object@target <- -0.3
  object@prob <- 1.1
  object@a <- -2
  object@b <- 0

  expect_equal(
    validate_stopping_lowest_dose_hsr_beta(object),
    c(
      "target must be probability > 0 and < 1",
      "prob must be probability > 0 and < 1",
      "Beta distribution shape parameter a must me a real number > 0",
      "Beta distribution shape parameter b must me a real number > 0"
    )
  )
})
