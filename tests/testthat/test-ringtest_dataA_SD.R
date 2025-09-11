
# test routine ---------------------------------------------------------------------
# Simulation with the GUTS model, compare against reference values.
# Dataset A is part of EFSA's *GUTS* ringtest to compare results from software
# implementations of [GUTS-RED models][GUTS-RED-models] (EFSA 2018). The
# ringtest focused on [GUTS-RED-IT][GUTS_RED_IT] and [GUTS-RED-SD][GUTS_RED_SD]
# models.
test_that("GUTS-RED-SD simulation", {
  tol <- 1e-4

  # GUTS SD
  guts_sd <- GUTS_RED_SD() %>%
    set_endpoints("S") %>%
    set_param(c(hb = 0, kd = 0.6, z = 5, kk = 0.5)) %>%
    set_noexposure() %>%
    set_bounds(list(
      kd = c(0.001, 1),
      z = c(0.01, 5),
      kk = c(0.01, 5),
      hb = c(0, 0.5)
    ))

  #####  optimise kD and TD parameters  ##########################
  # GUTS SD
  dt_control <- ringtest_sd %>%
    dplyr::filter(replicate == "Control") %>%
    dplyr::mutate(S = Nsurv / max(Nsurv)) %>%
    dplyr::select(time, S)

  fit <- calibrate(
    x = guts_sd,
    par = c(hb = 0),
    data = dt_control,
    output = "S",
    method = "Brent",
    verbose=FALSE
  )

  guts_sd <- guts_sd %>% set_param(c(hb = as.numeric(fit$par)))

  #####

  data <- ringtest_sd %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(S = Nsurv / max(Nsurv)) %>%
    dplyr::select(time, S, replicate, conc) %>%
    dplyr::ungroup()

  calib_par <- calibrate(guts_sd,
                         data = data,
                         par = c(kd = 1, z = 0.1, kk = 0.1),
                         output = "S",
                         maxsteps = 10^6,
                         approx = "linear",
                         hmax = 0.01,
                         verbose = FALSE,
                         hessian = TRUE
  )

  guts_sd <- guts_sd %>%
    set_param(calib_par$par) %>%
    set_noexposure()

  suppressMessages({
  res <- lik_profile(
    x = guts_sd,
    data = data,
    output = "S",
    par = calib_par$par,
    bounds = list(
      kd = c(0.001, 1),
      z = c(0.01, 10),
      kk = c(0.01, 5),
      hb = c(0, 1)
    ),
    refit = FALSE,
    type = "fine",
    method = "Brent",
    data_type = "continuous"
  )
  })

  # # tests for best fit parameters
  expect_equal(round(guts_sd@param$hb,6), 0.011794, tolerance = tol)
  expect_equal(round(guts_sd@param$kd,6), 0.685552, tolerance = tol)
  expect_equal(round(guts_sd@param$z,6), 3.003922, tolerance = tol)
  expect_equal(round(guts_sd@param$kk,6), 0.728990, tolerance = tol)
  expect_equal(res$kd$confidence_interval[1], 0.6671638, tolerance = tol)
  expect_equal(res$kd$confidence_interval[2], 0.7041692, tolerance = tol)
  expect_equal(res$kk$confidence_interval[1], 0.6889944, tolerance = tol)
  expect_equal(res$kk$confidence_interval[2], 0.7706515, tolerance = tol)
  expect_equal(res$z$confidence_interval[1], 2.964857, tolerance = tol)
  expect_equal(res$z$confidence_interval[2], 3.041375, tolerance = tol)
})
