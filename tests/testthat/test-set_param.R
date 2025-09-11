test_that("vanilla usage", {
  sc1 <- new("EffectScenario", name="foo")
  sc2 <- new("EffectScenario", name="foo", tag="bar")

  p <- list(kd=1, hb=2, alpha=3, beta=4)
  ps1 <- parameter_set("foo", p)

  ## vanilla scenarios
  # single scenario, set vector of atomic values
  expect_equal(set_param(sc1, unlist(p))@param, p)
  # single scenario, single parameter set
  expect_equal(set_param(sc1, ps1)@param, p)
  # multiple scenarios, atomic vector
  lst <- set_param(c(sc1,sc1), unlist(p))
  expect_equal(lst[[1]]@param, p)
  expect_equal(lst[[2]]@param, p)
  # multiple scenarios, parameter set
  lst <- set_param(c(sc1,sc1), ps1)
  expect_equal(lst[[1]]@param, p)
  expect_equal(lst[[2]]@param, p)
})

test_that("special cases", {
  sc1 <- minnow_it
  sc2 <- sc1 %>% set_times(sc1@times + max(sc1@times))
  p <- list(kd=1)
  ps1 <- parameter_set(sc1@name, p, tag=sc1@tag)

  # all scenarios within a sequence need to be modified
  sequence(seq=c(sc1, sc2)) %>%
    set_param(ps1) -> seq
  expect_equal(length(seq@scenarios), 2)
  expect_equal(seq@scenarios[[1]]@param$kd, p$kd)
  expect_equal(seq@scenarios[[2]]@param$kd, p$kd)
})

test_that("arg x=ScenarioSequence", {
  sc <- GUTS_RED_IT() %>% set_times(0:5)
  sq <- sequence(c(sc, sc), breaks=3)
  pv <- c(hb=23)
  ps <- parameter_set(sc@name, pv)

  ## single sequence
  # param=vector
  foo <- set_param(sq, pv)
  expect_equal(foo@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo@scenarios[[2]]@param, as.list(pv))
  # param=ParameterSet
  foo <- set_param(sq, ps)
  expect_equal(foo@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo@scenarios[[2]]@param, as.list(pv))

  ## multiple sequences
  lst <- list(sq, sq)
  # param=vector
  foo <- set_param(lst, pv)
  expect_equal(length(foo), 2)
  expect_equal(foo[[1]]@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo[[1]]@scenarios[[2]]@param, as.list(pv))
  expect_equal(foo[[2]]@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo[[2]]@scenarios[[2]]@param, as.list(pv))
  # param=ParameterSet
  foo <- set_param(lst, ps)
  expect_equal(length(foo), 2)
  expect_equal(foo[[1]]@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo[[1]]@scenarios[[2]]@param, as.list(pv))
  expect_equal(foo[[2]]@scenarios[[1]]@param, as.list(pv))
  expect_equal(foo[[2]]@scenarios[[2]]@param, as.list(pv))
})

test_that("arg x=CalibrationSet", {
  sc <- GUTS_RED_IT() %>% set_times(0:5)
  cs <- caliset(sc, data=data.frame(t=0, o=0))
  pv <- c(hb=23)
  ps <- parameter_set(sc@name, pv)

  ## single sequence
  # param=vector
  foo <- set_param(cs, pv)
  expect_equal(foo@scenario@param, as.list(pv))
  # param=ParameterSet
  foo <- set_param(cs, ps)
  expect_equal(foo@scenario@param, as.list(pv))

  ## multiple sequences
  lst <- list(cs, cs)
  # param=vector
  foo <- set_param(lst, pv)
  expect_equal(length(foo), 2)
  expect_equal(foo[[1]]@scenario@param, as.list(pv))
  expect_equal(foo[[2]]@scenario@param, as.list(pv))
  # param=ParameterSet
  foo <- set_param(lst, ps)
  expect_equal(length(foo), 2)
  expect_equal(foo[[1]]@scenario@param, as.list(pv))
  expect_equal(foo[[2]]@scenario@param, as.list(pv))
})

test_that("invalid arguments", {
  sc1 <- new("EffectScenario", name="foo")
  sc2 <- new("EffectScenario", name="foo", tag="bar")

  p <- list(kd=1)
  ps1 <- parameter_set("foo", p)
  ps2 <- parameter_set("foo", tag="bar", param=p)

  # multiple scenarios, one mismatch
  expect_error(set_param(c(sc1, sc2), ps1))
  # multiple parameter sets match, i.e. ambiguous assignments
  suppressMessages(expect_error(set_param(sc1, list(ps1, ps1))))
  # inconsistent types
  expect_error(set_param(sc1, list(1, ps1)))
  # model & tag dont match
  expect_error(set_param(sc1, ps2))
  # warn if invalid parameters were passed as argument
  sc3 <- new("EffectScenario", param.req=c("a"))
  expect_warning(set_param(sc3, c("b"=1)))
})
