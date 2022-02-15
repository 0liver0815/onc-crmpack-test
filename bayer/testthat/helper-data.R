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
