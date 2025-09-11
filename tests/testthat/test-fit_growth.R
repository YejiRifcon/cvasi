test_that("arg x=scenario (Lemna_Schmitt)", {
  sc <- metsulfuron %>% set_noexposure()
  df <- data.frame(t=0:2, o=12)
  td <- tox_data(df)
  # w/o arg data
  expect_error(fit_growth(sc), "data. is missing")

  # w/ arg data=data.frame
  fit_growth(sc, data=df)
  succeed()

  # w/ arg data=tox_data
  fit_growth(sc, data=td)
  succeed()

  # w/ arg data, but no exposure info anywhere. this can happen if no exposure
  # series was set for a newly created scenario.
  sc2 <- metsulfuron
  sc2@exposure@series <- data.frame(t=numeric(0), c=numeric(0))
  fit_growth(sc2, data=df)
  succeed()
})

test_that("arg x=scenario w/ exposure", {
  sc <- metsulfuron %>% set_exposure(data.frame(time=0:2, 1))
  df <- data.frame(t=0:2, o=12)

  # w/ arg data=data.frame
  expect_error(fit_growth(sc, data=df), "must not include exposure")
})

test_that("arg x=sequence", {
  sc <- metsulfuron %>% set_noexposure()
  sq <- sequence(list(sc, sc), breaks=1)
  df <- data.frame(t=0:2, o=12)

  # w/ arg data=data.frame
  fit_growth(sq, data=df)
  succeed()
})

test_that("arg x=caliset", {
  sc <- metsulfuron %>% set_noexposure()
  df <- data.frame(t=0:2, o=12)
  cs <- caliset(sc, df)

  expect_error(fit_growth(cs), "not supported")
})

test_that("arg x=list of calisets", {
  sc <- metsulfuron %>% set_noexposure()
  df <- data.frame(t=0:2, o=12)
  cs <- list(caliset(sc, df))

  fit_growth(cs)
  succeed()

  # w/ arg data=fails
  expect_error(fit_growth(cs, data=df), "not supported")
})

test_that("arg x=unsupported scenario", {
  sc <- new("EffectScenario") %>% set_noexposure()
  df <- data.frame(t=0:2, o=12)

  # scenario w/o data
  expect_error(fit_growth(sc), "data. is missing")
  # scenario w/ data
  expect_error(fit_growth(sc, data=df), "not supported yet")
  # caliset
  cs <- caliset(sc, data=df)
  expect_error(fit_growth(cs), "not supported")
})

# arg x=fit: does this need support?
