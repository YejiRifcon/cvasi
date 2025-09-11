test_that("arg sets", {
  # empty set
  rs <- eval_cs(list(), output="foo")
  expect_contains(names(rs), c("pred","obs","wgts","times","tags")) # partial check only
  expect_length(rs$obs, 0)
  expect_length(rs$pred, 0)
  expect_length(rs$times, 0)

  # single set
  sc <- GUTS_RED_IT()
  cs <- caliset(sc, data.frame(t=0:2, o=3:5), weight=6:8, tag="foo")
  pred <- 1:3
  rs <- with_mocked_bindings(eval_cs(list(cs), output="bar"),
    simulate=function(x, ...) {
      data.frame(time=get_times(x), "bar"=pred)
    })
  expect_equal(rs$obs, cs@data[, 2])
  expect_equal(rs$pred, pred)
  expect_equal(rs$times, cs@data[, 1])
  expect_equal(rs$wgts, cs@weight)
  expect_equal(rs$tags, as.list(rep(cs@tag, nrow(cs@data))))

  # multiple sets
  cs2 <- caliset(sc, data.frame(t=1:3, o=4:6), weight=7:9, tag="bar")
  rs <- with_mocked_bindings(eval_cs(list(cs, cs2), output="bar"),
    simulate=function(x, ...) {
     data.frame(time=get_times(x), "bar"=pred)
    })
  expect_equal(rs$obs, c(cs@data[, 2], cs2@data[, 2]))
  expect_equal(rs$pred, c(pred, pred))
  expect_equal(rs$times, c(cs@data[, 1], cs2@data[, 1]))
  expect_equal(rs$wgts, c(cs@weight, cs2@weight))
  expect_equal(rs$tags, as.list(c(rep(cs@tag, nrow(cs@data)),
                                rep(cs2@tag, nrow(cs2@data)))))
})

test_that("arg output", {
  sc <- GUTS_RED_IT()
  cs <- caliset(sc, data.frame(t=0:2, o=3))
  rs <- with_mocked_bindings(eval_cs(list(cs), output="baz"),
    simulate=function(x, ...) {
     data.frame(time=get_times(x), bar=8, baz=9)
    })
  expect_equal(rs$pred, c(9, 9, 9))
})

test_that("arg method/ode_method", {
  sc <- GUTS_RED_IT()
  cs <- caliset(sc, data.frame(t=0:1, o=3))
  fsim <- function(x, method=0, ...) data.frame(time=get_times(x), foo=method)

  # no method requested
  rs <- with_mocked_bindings(eval_cs(list(cs), output="foo"), simulate=fsim)
  expect_equal(rs$pred, c(0, 0))
  # arg 'method' set
  rs <- with_mocked_bindings(eval_cs(list(cs), output="foo", method="bar"), simulate=fsim)
  expect_equal(rs$pred, c("bar", "bar"))
  # arg 'ode_method' set
  rs <- with_mocked_bindings(eval_cs(list(cs), output="foo", ode_method="bar"), simulate=fsim)
  expect_equal(rs$pred, c("bar", "bar"))
  # arg 'ode_method' overwrites arg 'method'
  rs <- with_mocked_bindings(eval_cs(list(cs), output="foo", ode_method="bar", method="baz"), simulate=fsim)
  expect_equal(rs$pred, c("bar", "bar"))
})

test_that("arg verbose", {
  sc <- GUTS_RED_IT()
  cs <- caliset(sc, data.frame(t=0:1, o=3))
  fsim <- function(x, method=0, ...) data.frame(time=get_times(x), foo=NA)

  # message
  expect_message(with_mocked_bindings(eval_cs(list(cs), output="foo", verbose=TRUE), simulate=fsim), "contains NA")
  # no message
  expect_no_message(with_mocked_bindings(eval_cs(list(cs), output="foo", verbose=FALSE), simulate=fsim))
})

test_that("arg .suppress", {
  sc <- GUTS_RED_IT()
  cs <- caliset(sc, data.frame(t=0:1, o=3))
  ferror <- function(x, ...) stop("foobar")

  # error: try-error object
  expect_warning(
      with_mocked_bindings(rs <- eval_cs(list(cs), output="foo", verbose=FALSE, .suppress=FALSE), simulate=ferror),
      "with caliset #1"
  )
  expect_true(rs$is_error)
  expect_false(rs$is_issue)
  expect_equal(rs$err_msg, "foobar")

  # issue: desolve aborted
  faborted <- function(x, ...) {
    df <- data.frame(time=0:1, foo=0)
    attr(df, "desolve_diagn") <- list(istate=-1)
    df
  }
  with_mocked_bindings(rs <- eval_cs(list(cs), output="foo", verbose=FALSE, .suppress=TRUE), simulate=faborted)
  expect_false(rs$is_error)
  expect_true(rs$is_issue)
  expect_equal(rs$err_msg, "simulation terminated early")
})
