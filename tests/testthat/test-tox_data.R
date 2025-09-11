test_that("arg:data invalid", {
  expect_error(tox_data(list(1, 2)), "be a data.frame")
  expect_warning(tox_data(data.frame(t=0, n=1, t=2, c=3, foo=4)), "more than four columns")
  expect_error(tox_data(data.frame(t="foo", n=0)), "first column must be numeric")
  expect_error(tox_data(data.frame(t=0, n="foo")), "second column must be numeric")
  expect_error(tox_data(data.frame(t=0, n=1, t=2, conc="foo")), "fourth column must be numeric")
  expect_error(tox_data(data.frame(t=0, n=1, t=2, conc=3), exposure=list()), "Confusion about exposure")
})

test_that("arg:data", {
  df <- data.frame(t=0:4, n=1)

  # various types of data arguments
  # matrix
  dat <- tox_data(as.matrix(df))@data
  expect_equal(length(dat), 1) # list of trial data
  dat <- dat[[1]] # unnest list
  expect_true(is.data.frame(dat))
  expect_equal(length(dat), 2)
  expect_equal(names(dat), names(df))
  expect_equal(dplyr::select(dat, 1, 2), df, ignore_attr=TRUE)

  # data.frame
  dat <- tox_data(df)@data
  expect_equal(length(dat), 1) # list of trial data
  dat <- dat[[1]] # unnest list
  expect_true(is.data.frame(dat))
  expect_equal(length(dat), 2)
  expect_equal(names(dat), names(df))
  expect_equal(dplyr::select(dat, 1, 2), df, ignore_attr=TRUE)

  # tibble
  dat <- tox_data(tibble::as_tibble(df))@data
  expect_equal(length(dat), 1) # list of trial data
  dat <- dat[[1]] # unnest list
  expect_true(is.data.frame(dat))
  expect_equal(length(dat), 2)
  expect_equal(names(dat), names(df))
  expect_equal(dplyr::select(dat, 1, 2), df, ignore_attr=TRUE)
})

test_that("arg:data=data.frame format", {
  # basic data, two cols
  df2 <- data.frame(t=0:4, n=1)
  dat <- tox_data(df2)
  # observed
  expect_true(is(dat, "ToxData"))
  expect_true(is.list(dat@data))
  expect_equal(length(dat@data), 1)
  expect_equal(names(dat@data), "trial") # default trial id
  expect_true(is.data.frame(dat@data[[1]]))
  expect_equal(dat@data[[1]], df2, ignore_attr=TRUE)
  expect_equal(names(dat@data[[1]]), names(df2))
  # exposure
  expect_true(is.list(dat@exposure))
  expect_equal(length(dat@exposure), 0)
  # qty
  expect_equal(dat@qty_cols, "n")

  # with trial id, three cols
  df3 <- data.frame(t=0:4, n=1, t="foo")
  dat <- tox_data(df3)
  # single trial
  expect_equal(length(dat@data), 1)
  expect_equal(names(dat@data), "foo")
  expect_equal(dat@data[[1]], df3[, c(1, 2)], ignore_attr=TRUE)
  expect_equal(length(dat@exposure), 0)
  expect_equal(dat@qty_cols, "n")
  # multiple trials
  df3a <- data.frame(t=0:2, n=1, trial="foo")
  df3b <- data.frame(t=0:2, n=2, trial="bar")
  df3 <- rbind(df3a, df3b)
  dat <- tox_data(df3)

  expect_equal(length(dat@data), 2)
  expect_equal(names(dat@data), c("foo", "bar"))
  expect_equal(dat@data[["foo"]], df3a[, c(1, 2)], ignore_attr=TRUE)
  expect_equal(dat@data[["bar"]], df3b[, c(1, 2)], ignore_attr=TRUE)
  expect_false(tibble::is_tibble(dat@data[["foo"]]))
  expect_equal(length(dat@exposure), 0)
  expect_equal(dat@qty_cols, "n")

  # with concentration, four cols
  df4a <- data.frame(t=0:2, n=3, trial="foo", c=5)
  df4b <- data.frame(t=0:2, n=4, trial="bar", c=c(6, 6, 7))
  df4 <- rbind(df4a, df4b)
  dat <- tox_data(df4)

  expect_equal(length(dat@data), 2)
  expect_equal(names(dat@data), c("foo", "bar"))
  expect_equal(dat@data[["foo"]], df4a[, c(1, 2)], ignore_attr=TRUE)
  expect_equal(dat@data[["bar"]], df4b[, c(1, 2)], ignore_attr=TRUE)
  expect_false(tibble::is_tibble(dat@data[["foo"]]))
  expect_equal(length(dat@exposure), 2)
  expect_equal(names(dat@exposure), c("foo", "bar"))
  expect_equal(dat@exposure[["foo"]], df4a[, c(1, 4)], ignore_attr=TRUE)
  expect_equal(dat@exposure[["bar"]], df4b[, c(1, 4)], ignore_attr=TRUE)
  expect_equal(dat@qty_cols, "n")
})

test_that("arg:data=data.frame with NAs", {
  # basic data, two cols
  df2 <- data.frame(t=0:2, n=c(1, NA_real_, 2))
  dat <- tox_data(df2)
  # observed without NAs
  expect_equal(length(dat@data), 1)
  expect_equal(dat@data[[1]], df2 %>% dplyr::filter(!is.na(n)))

  # with concentration, four cols
  df4 <- data.frame(t=0:3, n=c(1, NA_real_, 2, 3), trial="foo", c=c(0, 1, NA_real_, 2))
  dat <- tox_data(df4)

  expect_equal(length(dat@data), 1)
  # observed without NAs
  expect_equal(dat@data[[1]], df4 %>% dplyr::filter(!is.na(n)) %>% dplyr::select(t, n))
  # exposure without NAs
  expect_equal(dat@exposure[[1]], df4 %>% dplyr::filter(!is.na(c)) %>% dplyr::select(t, c))
})

# TODO test tox_data with exposure argument
# TODO test surv_data

test_that("tox_exposure", {
  # named numeric vector/list
  vec <- c(foo=1)
  exp <- tox_exposure(vec)

  expect_equal(length(exp), 1)
  expect_equal(names(exp), names(vec))
  expect_equal(exp[[1]][, 2], vec[[1]])

  vec <- c(foo=1, bar=2)
  exp <- tox_exposure(vec)

  expect_equal(length(exp), 2)
  expect_equal(names(exp), names(vec))
  expect_equal(exp[[1]][, 2], vec[[1]])
  expect_equal(exp[[2]][, 2], vec[[2]])

  vec <- list(foo=1, bar=2)
  exp <- tox_exposure(vec)

  expect_equal(length(exp), 2)
  expect_equal(names(exp), names(vec))
  expect_equal(exp[[1]][, 2], vec[[1]])
  expect_equal(exp[[2]][, 2], vec[[2]])

  # data.frame w/o trial ids -> fails
  df <- data.frame(time=0, conc=1)
  expect_error(tox_exposure(df), "at least three columns")

  # data.frame with trial id
  df <- data.frame(time=0, conc=1, trial="foo")
  exp <- tox_exposure(df)

  expect_equal(length(exp), 1)
  expect_equal(names(exp), "foo")
  expect_equal(exp[[1]], df[, c(1, 2)], ignore_attr=TRUE)

  dfa <- data.frame(time=0:1, conc=1, trial="foo")
  dfb <- data.frame(time=0:1, conc=c(2, 3), trial="bar")
  df <- rbind(dfa, dfb)
  exp <- tox_exposure(df)

  expect_equal(length(exp), 2)
  expect_equal(sort(names(exp)), sort(unique(df$trial)))
  expect_equal(exp[["foo"]], dfa[, c(1, 2)], ignore_attr=TRUE)
  expect_equal(exp[["bar"]], dfb[, c(1, 2)], ignore_attr=TRUE)

  # list of various object types
  lst <- list(foo=1, bar=data.frame(t=1, c=2))
  exp <- tox_exposure(lst)

  expect_equal(length(exp), 2)
  expect_equal(sort(names(exp)), sort(unique(df$trial)))
  expect_equal(exp[["foo"]][, 2], lst$foo, ignore_attr=TRUE)
  expect_equal(exp[["bar"]], lst$bar[, c(1, 2)], ignore_attr=TRUE)
})

test_that("tox_exposure invalid arg", {
  expect_error(tox_exposure())
  # missing trial ids
  expect_error(tox_exposure(c(1, 2)), "empty or missing")
  expect_error(tox_exposure(list(1, 2)), "empty or missing")
  expect_error(tox_exposure(list(foo=1, 2)), "empty or missing")
  expect_error(tox_exposure(data.frame(time=0, conc=1, trial=c("foo","bar",""))), "empty or missing")

  # invalid data.frame format
  expect_error(tox_exposure(data.frame(t=0, c=0)), "at least three columns")
  expect_error(tox_exposure(data.frame(t="foo", c=0, tr="bar")), "must be numeric")
  expect_error(tox_exposure(data.frame(t=0, c="foo", tr="bar")), "must be numeric")

  # other invalid argument
  expect_error(tox_exposure(minnow_it), "invalid type")

  # invalid list elements
  expect_error(tox_exposure(list(foo=1, bar="a")), "invalid type")
})

test_that("td2cs", {
  ## valid arguments
  sc <- Lemna_SETAC() %>% set_times(0:4)
  # single tox_data dataset
  td <- tox_data(data.frame(time=0, obs=0, trial="foo"))
  cs <- td2cs(sc, td)
  expect_equal(length(cs), 1)
  expect_true(all(is_caliset(cs)))
  expect_equal(cs[[1]]@scenario, sc)
  expect_equal(cs[[1]]@data, td@data[[1]])
  expect_equal(cs[[1]]@tag, "foo")

  # scenario is a sequence, single tox_data dataset
  sq <- sequence(list(sc, sc), breaks=2)
  td <- tox_data(data.frame(time=0, obs=0, trial="foo"))
  cs <- td2cs(sq, td)
  expect_equal(length(cs), 1)
  expect_true(all(is_caliset(cs)))
  expect_equal(cs[[1]]@scenario, sq)
  expect_equal(cs[[1]]@data, td@data[[1]])
  expect_equal(cs[[1]]@tag, "foo")

  # multiple tox_data datasets
  td <- tox_data(data.frame(time=0, obs=0:1, trial=c("foo", "bar")))
  cs <- td2cs(sc, td)
  expect_equal(length(cs), 2)
  expect_true(all(is_caliset(cs)))
  expect_equal(cs[[1]]@scenario, sc)
  expect_equal(cs[[2]]@scenario, sc)
  expect_equal(cs[[1]]@data, td@data[[1]])
  expect_equal(cs[[2]]@data, td@data[[2]])
  expect_equal(cs[[1]]@tag, "foo")
  expect_equal(cs[[2]]@tag, "bar")

  # single tox_data dataset, multiple exposure series
  exp <- list("baz"=data.frame(t=0, c=0), "foo"=data.frame(t=0, c=1))
  td <- tox_data(data.frame(time=0, obs=0, trial="foo"), exposure=exp)
  cs <- td2cs(sc, td)
  expect_equal(length(cs), 1)
  expect_equal(cs[[1]]@scenario@param, sc@param) # use parameter list as a proxy, bc they do differ in exposure
  expect_equal(cs[[1]]@scenario@times, sc@times) # use parameter list as a proxy, bc they do differ in exposure
  expect_equal(cs[[1]]@data, td@data[[1]])
  expect_equal(cs[[1]]@tag, "foo")
  expect_equal(cs[[1]]@scenario@exposure@series, exp$foo)

  # single tox_data dataset, rename output var
  td <- tox_data(data.frame(time=0, foo=0))
  cs <- td2cs(sc, td, output_var="bar")
  expect_equal(length(cs), 1)
  expect_true(all(is_caliset(cs)))
  expect_equal(cs[[1]]@scenario, sc)
  expect_equal(cs[[1]]@data, data.frame(time=0, bar=0))

  ## invalid arguments
  sc <- Lemna_SETAC()
  expect_error(td2cs(), "scenario. is missing")
  expect_error(td2cs(list(sc)), "scenario. must be of length one")
  expect_error(td2cs("foo"), "scenario. must be a scenario")

  expect_error(td2cs(sc), "data. is missing")
  expect_error(td2cs(sc, data="foo"), "data. must be a tox_data object")

  expect_error(td2cs(sc, td, output_var=1), "must be a string")
  expect_error(td2cs(sc, td, output_var=c("foo", "bar")), "must be of length one")
})
