# maxDose-IncrementsHSRBeta ----

my_data <- h_get_data()
my_data@y[my_data@cohort==3] <- c(0L,0L,1L,1L)

test_that("IncrementsHSRBeta works correctly if toxcicity probability is below threshold probability", {
  increments <- IncrementsHSRBeta(target = 0.3, prob = 0.95)
  result <- maxDose(
    increments,
    data = my_data
  )
  expect_equal(result, 300) # maxDose is 300 as toxicity probability of no dose is above 0.95.
})

test_that("IncrementsHSRBeta works correctly if toxcicity probability is above threshold probability", {
  increments <- IncrementsHSRBeta(target = 0.3, prob = 0.9)
  result <- maxDose(
    increments,
    data = my_data
  )
  expect_equal(result, 75) # maxDose is 75 as toxicity probability of dose 100 is above 0.90.
})


my_data@y[my_data@cohort==1] <- c(0L,1L,1L,1L)

test_that("IncrementsHSRBeta works correctly if toxcicity probability of first active dose is above threshold probability", {
  increments <- IncrementsHSRBeta(target = 0.3, prob = 0.95)
  result <- maxDose(
    increments,
    data = my_data
  )
  expect_equal(result, 25) # maxDose is 25 as toxicity probability of dose 25 is above 0.95 and placebo used.
})


my_data <- h_get_data()
my_data@y[my_data@x==0.001] <- c(1L,1L,1L)

test_that("IncrementsHSRBeta works correctly if toxcicity probability of placebo is above threshold probability", {
  increments <- IncrementsHSRBeta(target = 0.3, prob = 0.95)
  result <- maxDose(
    increments,
    data = my_data
  )
  expect_equal(result, 300) # maxDose is 300 as placebo is ignored.
})
