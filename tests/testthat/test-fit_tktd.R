test_that("dynamic call routing: list of calisets", {
  # create dummy caliset
  cs <- caliset(scenario=Lemna_SETAC(), data=data.frame(time=0:1, obs=0))
  my_tktd <- function(x, ...) { return(TRUE) }

  # check if dynamic call is routed correctly
  rs <- FALSE
  with_mocked_bindings(
    rs <- fit_tktd(list(cs)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(rs)
})

test_that("dynamic call routing: scenario sequence", {
  # create dummy sequence
  sc <- Lemna_SETAC() %>% set_times(0:8)
  sq <- sequence(list(sc, sc), breaks=3)
  my_tktd <- function(x, ...) { return(TRUE) }

  # check if dynamic call is routed correctly
  rs <- FALSE
  with_mocked_bindings(
    rs <- fit_tktd(sq, data=data.frame(time=0:1, obs=0)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(rs)
})

test_that("arg x=scenario", {
  ## Valid arguments
  # scenario with data
  sc <- Lemna_SETAC() %>% set_times(0:8)
  my_tktd <- function(x, ...) { return(x) }

  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(sc, data=data.frame(time=0:1, obs=0)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(is.list(rs))
  expect_true(is_caliset(rs[[1]]))
  expect_equal(rs[[1]]@scenario, sc)

  # sequence with data
  sq <- sequence(list(sc, sc), breaks=3)

  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(sq, data=data.frame(time=0:1, obs=0)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(is.list(rs))
  expect_true(is_caliset(rs[[1]]))
  expect_equal(rs[[1]]@scenario, sq)

  ## Invalid arguments
  # list of scenarios with data
  expect_error(fit_tktd(list(sc), data=data.frame(time=0:1, obs=0)), "not supported")
  # scenario without data
  expect_error(fit_tktd(sc), "is a scenario")

})

test_that("arg x=cvasi_fit", {
  fit <- list(par=c("k_photo_max"=0.123456))
  class(fit) <- c("cvasi_fit", "list")
  cs <- caliset(Lemna_SETAC(), data.frame(time=0:1, obs=0))
  csref <- caliset(Lemna_SETAC() %>% set_param(fit$par), data.frame(time=0:1, obs=0))
  my_tktd <- function(x, ...) { return(x) }

  ## Valid arguments
  # data=single caliset
  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(fit, data=cs),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_equal(rs, list(csref))

  # data=list of calisets
  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(fit, data=list(cs, cs)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_equal(rs, list(csref, csref))

  ## Invalid arguments
  # data=missing
  expect_error(fit_tktd(fit), "data. is missing")
  # data=data.frame
  expect_error(fit_tktd(fit, data=data.frame(time=0:2, obs=0)), "only contain caliset")
})

test_that("arg x=caliset(s)", {
  ## Valid arguments
  sc <- Lemna_SETAC() %>% set_times(0:8)
  cs <- caliset(sc, data.frame(time=0:1, obs=0))
  my_tktd <- function(x, ...) { return(x) }

  # single caliset
  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(cs),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(is.list(rs))
  expect_equal(length(rs), 1)
  expect_true(is_caliset(rs[[1]]))
  expect_equal(rs[[1]]@scenario, sc)

  # multiple calisets
  rs <- NULL
  with_mocked_bindings(
    rs <- fit_tktd(list(cs, cs)),
    fit_tktd_lemna_setac = my_tktd
  )
  expect_true(is.list(rs))
  expect_equal(length(rs), 2)
  expect_equal(rs, list(cs, cs))

  ## Invalid arguments
  # caliset with data
  expect_error(fit_tktd(cs, data=data.frame(t=0, o=0)), "alisets.*cannot be used together")
})

test_that("arg x=unsupported model", {
  source(test_path("dummy.R"), local = TRUE)
  sc <- new("DummyScenario")
  df <- data.frame(t=0, obs=0)

  # args: scenario + data
  expect_error(fit_tktd(sc, data=df), "Scenario type not supported")

  # args: caliset
  cs <- caliset(sc, data=df)
  expect_error(fit_tktd(cs), "Scenario type not supported")
})
