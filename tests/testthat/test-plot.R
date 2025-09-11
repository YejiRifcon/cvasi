test_that("cvasi_drc", {
  df <- dose_response(minnow_it)
  plot(df)
  succeed()
})

test_that("cvasi_simulate", {
  df <- simulate(minnow_it)
  plot(df)
  succeed()
})

test_that("lik_profile", {
  data <- schmitt2013 %>% dplyr::filter(trial == "T0")

  # update metsulfuron
  scenario <- metsulfuron %>%
    set_param(c(k_phot_fix = TRUE, Emax = 1)) %>%
    set_init(c(BM = 12)) %>%
    set_bounds(list(k_phot_max=c(0, 1)))

  # Likelihood profiling
  suppressMessages(
    res <- lik_profile(
    x = td2cs(scenario, tox_data(data)),
    data = data,
    output = "BM",
    par = c("k_phot_max"=0.43),
    bounds = list(
      k_phot_max = list(0, 1)
    ),
    refit=FALSE
  ))
  # plot
  plot(res)
  succeed()
})
