# stopTrial-StoppingLowestDoseHSRBeta ----

# Sample data to test Stopping Rule lowest active dose is toxic.
my_data <- h_get_data()
my_model <- LogisticKadane(0.3, xmin = 0.001, xmax = 100)
my_options <- McmcOptions(
  burnin = 1, step = 1, samples = 1, rng_kind = "Mersenne-Twister", rng_seed = 94
)
my_samples <- mcmc(my_data, my_model, my_options)

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not toxic", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.9)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    FALSE,
    message = "Probability that the lowest active dose of 25 being toxic based on posterior Beta distribution using a Beta(1,1) prior is 24% and thus below the required 90% threshold."
  )
  expect_identical(result, expected) # Prob being toxic is 24% < 90%.
})

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not toxic", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.1)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    TRUE,
    message = "Probability that the lowest active dose of 25 being toxic based on posterior Beta distribution using a Beta(1,1) prior is 24% and thus above the required 10% threshold."
  )
  expect_identical(result, expected) # Prob being toxic is 24% > 10%.
})


h_get_data_no_plcb <- function() {
  x <- c(25, 25, 25, 50, 50, 50, 100, 100, 100)
  dose_grid <- c(seq(25, 300, 25))

  Data(
    x = x,
    y = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L),
    doseGrid = dose_grid,
    placebo = FALSE,
    ID = 1:9,
    cohort = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  )
}

my_data <- h_get_data_no_plcb()
my_model <- LogisticKadane(0.3, xmin = 0.001, xmax = 100)
my_options <- McmcOptions(
  burnin = 1, step = 1, samples = 1, rng_kind = "Mersenne-Twister", rng_seed = 94
)
my_samples <- mcmc(my_data, my_model, my_options)

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not toxic", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.9)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    FALSE,
    message = "Probability that the lowest active dose of 25 being toxic based on posterior Beta distribution using a Beta(1,1) prior is 24% and thus below the required 90% threshold."
  )
  expect_identical(result, expected) # Prob being toxic is 24% < 90%.
})

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not toxic", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.1)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    TRUE,
    message = "Probability that the lowest active dose of 25 being toxic based on posterior Beta distribution using a Beta(1,1) prior is 24% and thus above the required 10% threshold."
  )
  expect_identical(result, expected) # Prob being toxic is 24% > 10%.
})

my_data <- h_get_data()
my_data@x[my_data@cohort==1] <- c(0.001,75,75,75)

my_model <- LogisticKadane(0.3, xmin = 0.001, xmax = 100)
my_options <- McmcOptions(
  burnin = 1, step = 1, samples = 1, rng_kind = "Mersenne-Twister", rng_seed = 94
)
my_samples <- mcmc(my_data, my_model, my_options)

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not applied", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.1)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    FALSE,
    message = "Lowest active dose not tested, stopping rule not applied."
  )
  expect_identical(result, expected) # First active dose not applied.
})


my_data <- h_get_data_no_plcb()
my_data@x[my_data@cohort==1] <- c(75,75,75)

my_model <- LogisticKadane(0.3, xmin = 0.001, xmax = 100)
my_options <- McmcOptions(
  burnin = 1, step = 1, samples = 1, rng_kind = "Mersenne-Twister", rng_seed = 94
)
my_samples <- mcmc(my_data, my_model, my_options)

test_that("StoppingLowestDoseHSRBeta works correctly if first active dose is not applied", {
  stopping <- StoppingLowestDoseHSRBeta(target = 0.3, prob = 0.1)
  result <- stopTrial(
    stopping = stopping,
    dose = 300,
    samples = my_samples,
    model = my_model,
    data = my_data
  )
  expected <- structure(
    FALSE,
    message = "Lowest active dose not tested, stopping rule not applied."
  )
  expect_identical(result, expected) # First active dose not applied.
})


