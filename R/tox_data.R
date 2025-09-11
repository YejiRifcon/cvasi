# Class definition
#' @aliases tox_data
setClass("ToxData",
  slots = list(data="ANY",
               exposure="ANY",
               qty_cols="character"),
  prototype = list(data=list(),
                   exposure=list(),
                   qty_cols=character(0))
)

# TODO overload plot() to display tox data
# TODO overload show to have some basic stats
# TODO overload as.data.frame() to return the tox data (we ignore exposure for now)
# TODO examples
#
#' Prepare ecotox study data for fitting
#'
#' Takes ecotox study data in *long-form* tabular format and prepares it for
#' parameter fitting. It supports extracting (optional) exposure concentration from
#' tabular data (useful for e.g. studies of acute toxicity) or exposure series
#' can also provided as individual time-series.
#'
#' # Tabular format
#' The long-form tabular data must have at least two and at most four columns.
#' The position of the columns define what they represent, the column names
#' are ignored:
#'
#' - First column: time
#' - Second column: observed quantity, e.g. number of individuals
#' - (optional) Third column: Trial or treatment ID
#' - (optional) Fourth column: Concentration
#'
#' The first two columns, time and observed quantity, must always be present.
#' The third column, trial ID, is used to split the table by treatment so
#' that trials can later be handled individually. The fourth column, concentration,
#' can be used to also define the exposure level during the experiments.
#'
#' # Explicit exposure series
#' As an alternative to defining concentrations along observed data, exposure
#' can also be passed as a list of exposure levels and series with
#' argument `exposure`. It can be used to provide exposure series for each
#' trial. The following object types are supported to define exposure:
#'
#' - Numerical constants
#' - Tabular data, e.g. `data.frame`s
#' - [exposure series][ExposureSeries] objects
#'
#' If exposure is constant over time, exposure can be defined using a single
#' constant value. More complex exposure time-series can be defined using e.g.
#' `data.frame`s. Tabular data must have two columns with the first column
#' representing time and the second column representing exposure/concentrations.
#'
#'
#' @param data a `data.frame` with at least two and at most four columns; the first
#'   column must represent time, the second an observed quantity, the optional third
#'   a trial or treatment ID, and the optional fourth the concentration during
#'   the experiment
#' @param exposure an optional named list; names must correspond to trial IDs used for
#'   the `data` argument; values can be numeric constants, data.frames, or an
#'   [ExposureSeries] object
#' @return a `ToxData` object
#' @aliases ToxData-class
#' @export
#' @examples
#' library(dplyr)
#'
#' mydata <- schmitt2013 %>% tox_data()
tox_data <- function(data, exposure=NULL) {
  if(missing(data))
    stop("Argument 'data' is missing")
  if(is.matrix(data))
    data <- as.data.frame(data)
  if(!is.data.frame(data))
    stop("Argument 'data' must be a data.frame")

  # additional checks on `data`
  errs <- c()
  has_conc <- ncol(data) >= 4
  has_trial <- ncol(data) >= 3
  if(ncol(data) > 4)
    warning("Argument `data` has more than four columns, additional columns will be ignored")
  if(ncol(data) == 0) {
    errs <- c(errs, "Argument 'data' is empty")
  }
  if(ncol(data) >= 1) {
    if(!is.numeric(dplyr::pull(data, 1))) # time
      errs <- c(errs, "Argument 'data': first column must be numeric")
  }
  if(ncol(data) >= 2) {
    if(!is.numeric(dplyr::pull(data, 2))) # quantity
      errs <- c(errs, "Argument 'data': second column must be numeric")
  }
  # no check of trial id, yet
  if(has_conc) {
    if(!is.numeric(dplyr::pull(data, 4))) # concentrations
      errs <- c(errs, "Argument 'data': fourth column must be numeric")

    # if `exposure` argument is present AND `data` contains concentrations,
    # then it's unclear what to do
    if(!is.null(exposure)) {
      stop("Confusion about exposure series: argument 'data' contains concentrations, but argument 'exposure' is also set")
    }
    # select exposure data as time, conc, trial id
    exposure <- dplyr::select(data, 1, 4, 3)
    # drop concentration from observed data
    data <- dplyr::select(data, 1, 2, 3)
  }
  if(!has_trial) {
    # assign default trial name
    data <- dplyr::mutate(data, trial="trial")
    has_trial <- TRUE
  }

  if(length(errs) > 0) {
    stop("\n* ", paste(errs, collapse="\n* "))
  }

  # convert inputs to named list of exposure series
  if(!is.null(exposure)) {
    exposure <- tox_exposure(exposure)
  } else {
    exposure <- list()
  }

  nms <- names(data)
  # split into list of observations per trial
  if(has_trial) {
    data <- data %>%
      dplyr::mutate(dplyr::across(3, ~ factor(., levels=unique(.)))) %>%
      dplyr::filter(dplyr::if_any(2, ~ !is.na(.))) %>%
      dplyr::group_by(dplyr::pick(3))
    ids <- dplyr::group_keys(data) %>% dplyr::pull(1)
    # sort by time and split by trial id
    data <- data %>%
      dplyr::arrange(dplyr::pick(1)) %>%
      dplyr::group_split(.keep=FALSE)
    data <- setNames(as.list(data), ids)
    for(i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]]) # get rid of tibbles
    }
  }
  # table contains data of a single trial
  else {
    data <- list(dplyr::arrange(data, dplyr::pick(1))) # order by time
  }

  new("ToxData", data=data, exposure=exposure, qty_cols=nms[2])
}

# Function for user convenience, intended to work similar to [morse::survData()]
#' @autoglobal
surv_data <- function(data) {
  # coerce argument to data.frame, e.g. to convert morse's data types such as `survDataCstExp`
  data <- as.data.frame(data)

  # TODO check that columns exist
  # TODO check data types

  # reorder and select columns
  qty <- data %>% dplyr::select(time, Nsurv, replicate)
  expo <- data %>% dplyr::select(time, conc, replicate)

  tox_data(qty, expo)
}

# Converts argument to a named list of objects that can be used as
# exposure series, i.e. as arguments to [set_exposure()]. Element names
# will be used as trial ids and vice versa in order to later associate
# exposure series to observed data sets.
#
# @param exp a numeric vector, a data.frame, an exposure series, or a (named) list
#   of any of the aforementioned types
# @return named list of [ExposureSeries] objects
tox_exposure <- function(exp) {
  if(missing(exp))
    stop("Argument 'exp' is missing")
  # argument must either be a named vector/list or a data.frame with three columns
  if(is.vector(exp)) {
    # check that names exist and none are empty
    nms <- names(exp)
    if(is.null(nms) | any(nms == ""))
      stop("Some trial ids are empty or missing")
    exp <- as.list(exp)
  }
  else if(is.data.frame(exp)) {
    exp <- as.data.frame(exp) # get rid of tibbles
    if(ncol(exp) < 3)
      stop("Argument must have at least three columns")
    if(!is.numeric(exp[, 1]))
      stop("First column must be numeric")
    if(!is.numeric(exp[, 2]))
      stop("Second column must be numeric")

    # split by trial id
    exp <- exp %>%
      dplyr::mutate(dplyr::across(3, ~ factor(., levels=unique(.)))) %>%
      dplyr::filter(dplyr::if_any(2, ~ !is.na(.))) %>%
      dplyr::group_by(dplyr::pick(3))
    lst <- dplyr::group_split(exp, .keep=FALSE)
    nms <- dplyr::group_keys(exp) %>% dplyr::pull(1)
    exp <- setNames(as.list(lst), nms) # get rid of tibble prototype
    for(i in seq_along(exp)) {
      exp[[i]] <- as.data.frame(exp[[i]]) # get rid of tibble format
    }
  }
  else {
    stop("Argument is of invalid type")
  }

  # check if list elements have valid types, convert where possible
  for(i in seq_along(exp)) {
    elem <- exp[[i]]
    nm <- names(exp)[[i]]
    # constant value
    if(is.numeric(elem)) {
      if(length(elem) != 1)
        stop("Exposure value of trial '", nm, "' must have length one")
      exp[[i]] <- data.frame(time=0, conc=elem)
    }
    # data.frame
    else if(is.data.frame(elem)) {
      if(ncol(elem) != 2)
        stop("Exposure series of trial '", nm, "' must have two columns")
    }
    # exposure series
    else if(all(is_exp_series(elem)) & length(elem) == 1) {
      # no checks for now
    } else {
      stop("Exposure series of trial '", nm, "' has invalid type")
    }
  }

  if(length(exp) == 0)
    stop("Argument does not contain any exposure information")

  # check names
  nms <- names(exp)
  if(is.null(nms)) {
    nms <- c("")
  }
  if(any(nms == "")) {
    stop("Some trial ids are empty or missing")
  }
  if(any(table(nms) > 1)) {
    stop("Trial ids are not unique")
  }

  exp
}

#' Create *calisets* from *tox_data*
#'
#' This is a convenience function which eases the creation of [caliset]s
#' from a base scenario and a [tox_data] object. The scenario
#' will be used as is. If exposure series are defined by the `tox_data`
#' object, then these will be assigned to the scenario(s) accordingly.
#'
#' @param scenario a base [scenario]
#' @param data return value of [tox_data()]
#' @param output_var optional *character*, rename observed data column to this value
#' @return list of [caliset]s
#' @export
#' @examples
#' # Import trial data from Schmitt et al. (2013), including exposure
#' mydata <- tox_data(schmitt2013)
#'
#' # Example trial contained in data set: trial 'T0.32'
#' mydata@data[["T0.32"]]        # observed quantities
#' mydata@exposure[["T0.32"]]    # associated exposure series
#'
#' # Create list of calisets from full data set
#' lst <- td2cs(Lemna_SETAC(), mydata)
#'
#' # Example caliset representing conditions of trial 'T0.32'
#' lst[[2]]
td2cs <- function(scenario, data, output_var=NULL) {
  if(is.list(scenario))
    stop("Argument 'scenario' must be of length one")
  if(!is_scenario(scenario) & !is_sequence(scenario))
    stop("Argument 'scenario' must be a scenario")

  if(missing(data))
    stop("Argument 'data' is missing")
  if(!is(data, "ToxData"))
    stop("Argument 'data' must be a tox_data object")

  if(!is.null(output_var)) {
    if(!is.character(output_var))
      stop("Argument 'output_var' must be a string")
    if(length(output_var) != 1)
      stop("Argument 'output_var' must be of length one")
  }

  lst <- list()
  nms <- names(data@data)
  for(i in seq_along(data@data)) {
    # create a copy of the scenario, as it might get modified
    sc <- scenario
    nm <- nms[[i]]

    # make sure data series is a `data.frame` and rename column with observed data
    df <- as.data.frame(data@data[[i]])
    # if the name of the observed variable is defined
    if(!is.null(output_var)) {
      df <- setNames(df, c(names(df)[[1]], output_var))
    }
    # look if an exposure series exists for that trial id
    if(nm %in% names(data@exposure)) {
      sc <- set_exposure(sc, data@exposure[[nm]], reset_times=FALSE)
    }
    lst[[i]] <- caliset(sc, data=df, tag=nm)
  }
  lst
}
