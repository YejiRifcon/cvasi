
# test routine ---------------------------------------------------------------------
# Simulation with the GUTS model, compare against reference values.
# Dataset A is part of EFSA's *GUTS* ringtest to compare results from software
# implementations of [GUTS-RED models][GUTS-RED-models] (EFSA 2018). The
# ringtest focused on [GUTS-RED-IT][GUTS_RED_IT] and [GUTS-RED-SD][GUTS_RED_SD]
# models.
test_that("GUTS-RED-IT simulation", {
  tol <- 1e-4

  # define GUTS model
  guts_it <- GUTS_RED_IT() %>%
    set_endpoints("S") %>%
    set_param(c(hb = 0, kd = 0.5, alpha = 5, beta = 5)) %>%
    set_noexposure() %>%
    set_bounds(list(
      kd = list(0.001, 1),
      alpha = list(0.1, 10),
      beta = list(0.1, 10),
      hb = list(0, 1)
    ))

  #####  optimise kD and TD parameters  ##########################

  # fit background mortality hb
  # GUTS IT
  dt_control <- ringtest_it %>%
    dplyr::filter(replicate == "Control") %>%
    dplyr::mutate(S = Nsurv / max(Nsurv)) %>%
    dplyr::select(time, S)

  fit <- calibrate(
    x = guts_it,
    par = c(hb = 0),
    data = dt_control,
    output = "S",
    method = "Brent",
    verbose=FALSE
  )

  guts_it <- guts_it %>% set_param(c(hb = as.numeric(fit$par)))

  #####

  data <- ringtest_it %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(S = Nsurv / max(Nsurv)) %>%
    dplyr::select(time, S, replicate, conc) %>%
    dplyr::ungroup()

  calib_par <- calibrate(guts_it,
                         data = data,
                         par = c(kd = 0.5, alpha = 5, beta = 5),
                         output = "S",
                         maxsteps = 10^6,
                         approx = "linear",
                         hmax = 0.01,
                         verbose = FALSE,
                         hessian = TRUE
  )

  guts_it <- guts_it %>%
    set_param(calib_par$par)

  suppressMessages({
  res <- lik_profile(
    x = guts_it,
    data = data,
    par = calib_par$par,
    output = "S",
    pars_bound <- list(
      hb = list(0, 1),
      kd = list(0.001, 1),
      alpha = list(0.1, 5),
      beta = list(0.1, 5)
    ),
    hmax = 0.01,
    maxsteps = 10^6,
    data_type = "continuous"
  )
  })

  # # tests for best fit parameters
  expect_equal(round(guts_it@param$hb,6), 0.044089, tolerance = tol)
  expect_equal(round(guts_it@param$kd,6), 0.774101, tolerance = tol)
  expect_equal(round(guts_it@param$alpha,6), 5.563681, tolerance = tol)
  expect_equal(round(guts_it@param$beta,6), 5.085245, tolerance = tol)
  expect_equal(res$kd$confidence_interval[1], 0.7076738, tolerance = tol)
  expect_equal(res$kd$confidence_interval[2], 0.8490582, tolerance = tol)
  expect_equal(res$alpha$confidence_interval[1], 5.451594, tolerance = tol)
  expect_equal(res$alpha$confidence_interval[2], 5.563681, tolerance = tol)
  expect_equal(res$beta$confidence_interval[1], 4.514981, tolerance = tol)
  expect_equal(res$beta$confidence_interval[2], 5.085245, tolerance = tol)
})
