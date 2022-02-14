# StoppingLowestDoseHSRBeta-class ----

test_that(".StoppingLowestDoseHSRBeta works as expected", {
  result <- expect_silent(.StoppingLowestDoseHSRBeta())
  expect_valid(result, "StoppingLowestDoseHSRBeta")
})

# StoppingLowestDoseHSRBeta-constructor ----

test_that("StoppingLowestDoseHSRBeta object can be created with user constructor", {
  result <- expect_silent(StoppingLowestDoseHSRBeta())
  expect_valid(result, "StoppingLowestDoseHSRBeta")
})
