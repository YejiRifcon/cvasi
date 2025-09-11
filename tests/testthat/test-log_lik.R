test_that("log-likelihood works: continous data, incl. log scale", {

  npars <- 2
  obs <- c(12, 38, 92, 176, 176, 627, 1283, 2640)
  pred <- c(5.0, 21.256, 55.586, 144.162, 144.845, 574.043, 1323.999, 2632.258)
  res <- log_lik(npars = npars,
                 obs = obs,
                 pred = pred)
  res_log <- log_lik(npars = npars,
                   obs = obs,
                   pred = pred,
                   log_scale = TRUE)

  expect_equal(res, -39.074, tolerance = 10e-3)  # expected value comes from previous cvasi runs
  expect_equal(res_log, -4.50445, tolerance = 10e-3) # expected value comes from previous cvasi runs

})

test_that("log-likelihood works: count data", {

  # observational data
  dt <- ringtest_c %>%
    dplyr::filter(replicate == "E")
  myexposure <- dt %>%
    dplyr::select(time, conc)
  obs <- dt %>%
    dplyr::mutate(S=Nsurv / max(Nsurv)) %>%
    dplyr::select(time, S)
  # GUTS model
  GUTS_RED_IT() %>%
    set_param(c(hb=0)) %>%
    set_exposure(myexposure) -> myscenario
  # fit
  suppressWarnings({
    fit <- calibrate(
      x = myscenario,
      par = c(kd=1.2, alpha=9.2, beta=4.3),
      data = obs,
      output = "S",
      verbose=FALSE)
  })
  # update
  myscenario <- myscenario %>%
    set_param(fit$par)
  # simulate
  pred <- myscenario %>%
    simulate()
  # calc likelihood
  res <- log_lik(obs = obs$S,
          pred =  pred$S,
          data_type = "count")

  expect_equal(res, -0.21, tolerance = 10e-3) # expected value comes from previous cvasi runs
})



test_that("log-likelihood gives expected error", {

  npars <- 2
  obs <- c(12, 38, 92, 176, 176, 627, 1283, 2640)
  pred <- c(5.0, 21.256, 55.586, 144.162, 144.845, 574.043, 1323.999, 2632.258)

  expect_error(log_lik(npars = npars,
                       obs = obs[-1],
                       pred = pred))

  expect_error(log_lik(npars = "test",
                       obs = obs,
                       pred = pred))

  expect_error(log_lik(npars = npars,
                       obs = obs,
                       pred = pred,
                       data_type = "foo"))

})
