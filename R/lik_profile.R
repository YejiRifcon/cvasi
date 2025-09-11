# TODO rename and overload profile() function
# TODO remove parameter 'bounds'
# TODO accept a cvasi_fit object as arg 'x'

#' Likelihood profiling
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The aim of the function is two-fold: 1) estimate a 95% confidence
#' around each parameter of a calibrated model, and 2) see if perhaps a local minimum was found rather
#' than a global minimum. To achieve this, the likelihood profiling goes through
#' every parameter one by one. For each parameter,
#' the model is sequentially refit with the parameter value set to
#' increasingly lower and higher values, and the likelihood of the model given the
#' data calculated (using [log_lik()]). The likelihood is then compared
#' to the likelihood of the original model (using a likelihood ratio). This leads
#' to the development of a likelihood profile, from which a plot a 95%
#' confidence interval for the parameter is derived.
#'
#' The idea of the function is a variable stepwise algorithm:
#' When the likelihood ratio changes very little (less than `l_crit_min`), the stepsize is
#' increased (up to a maximum, specified by `f_step_max`). When the lik.
#' ratio changes too much (more than `l_crit_max`), the algorithm tries again
#' with a smaller stepsize (also bound to a minimum: `f_step_min`). Note that
#' the stepsize is used as a fraction of the parameter value that is tried.
#' To prevent very small stepsizes when the value goes towards zero (as can
#' be the case for effect thresholds), an absolute minimum
#' stepsize (`f_step_abs`), which is specified as a fraction of the best
#' parameter value (`Xhat`) (unless it is zero, then the algorithm takes
#' something small).
#'
#' Note that the likelihood of the model given the data can be calculated across
#' all datasets provided in the calibration set `x`, or calculated separately for
#' each individual dataset before being combined into one likelihood (by adjusting
#' the optional parameter `individual`). The latter
#' has the advantage that different datasets can be given different weights in
#' the likelihood calculation (using the "weight" slot of the [caliset] objects,`x`).
#' Further, for continuous data (e.g. biomass), the likelihood considers the variance (standard
#' deviation) in the log likelihood calculation, which can vary between datasets
#' when the likelihood is calculated for each dataset separately before combining
#' into an overall likelihood. The latter could be relevant when factors might lead
#' to variability between datasets (e.g. different labs, different animal culture,...)
#'
#' To conduct the likelihood calculations on separate datasets, the parameter `individual`
#' which by default is 'FALSE' can be set to 'TRUE'. Then, then log likelihoods
#' are calculated for each dataset individually (or in subgroups, using the "tag" names of the
#' [caliset] object, if provided, to group datasets with the same "tag"
#' before calculating the log likelihood). Subsequently, the log
#' likelihoods for the subsets are combined into an overall likelihood (considering
#' the *set* weights provided in the "weight" slot of the [caliset] object).
#' Note that for each *set* only 1 weight can be provided (i.e. not individual
#' weights for each datapoint within the *set*), and that *set* with the same tag
#' should have identical weight.
#'
#' The function was inspired by a MatLab BYOM v.6.8 procedure, created by
#' Tjalling Jager. For details, please refer to BYOM (http://debtox.info/byom.html)
#' as well as Jager (2021).
#'
#' @references
#' Jager T, 2021: Robust Likelihood-Based Optimization and Uncertainty Analysis
#' of Toxicokinetic-Toxicodynamic Models. Integrated Environmental Assessment and
#' Management 17:388-397. \doi{10.1002/ieam.4333}
#'
#' @param x either a single [scenario] or a list of [caliset] objects
#' @param par named vector - parameters (names and values) to be profiled
#' @param output character vector, name of output column of [simulate()] that
#'  is used in calibration
#' @param data only needed if `x` is a [scenario]
#' @param bounds optional list of lists (including lower and upper bound): uses defaults in `x` object, but
#'  can be overwritten here (e.g. bounds <- list(k_resp = list(0,10), k_phot_max = list(0,30)) )
#' @param refit if `TRUE` (default), refit if a better minimum is found
#' @param type `"fine"` or `"coarse"` (default) likelihood profiling
#' @param individual if `FALSE` (default), the log likelihood is calculated across
#' the whole dataset. Alternatively, if `TRUE`, log likelihoods are calculated for
#' each (group of) *set*(s) individually.
#' @param break_prof If `TRUE`, then stop the profiling if a better optimum is located.
#'    Default is `FALSE`.
#' @param log_scale `FALSE` (default), option to calculate the log likelihood on a
#' log scale (i.e., observations and predictions are log transformed during calculation)
#' @param data_type Character argument, `"continuous"` (default) or `"count"`, to specify the data type
#' for the log likelihood calculations.
#' @param ... additional parameters passed on to [calibrate()] and [simulate()]. To avoid
#'   parameter confusion, use argument `method` to select optimization algorithms
#'   of `calibrate()` and argument `ode_method` to select numerical integration
#'   schemes of package `deSolve`.
#' @autoglobal
#' @examples
#' # Example with Lemna model - physiological params
#' library(dplyr)
#'
#' # observations - control run
#' obs <- schmitt2013 %>%
#'   filter(trial == "T0")
#'
#' # update metsulfuron
#' myscenario <- metsulfuron %>%
#'   set_param(c(k_phot_fix = TRUE, Emax = 1)) %>%
#'   set_init(c(BM = 0.0012)) %>%
#'   set_noexposure() %>%
#'   set_bounds(list(k_phot_max=c(0, 1)))
#'
#' fit <- calibrate(
#'   x = myscenario,
#'   par = c(k_phot_max = 1),
#'   data = obs,
#'   output = "FrondNo",
#'   method = "Brent"
#' )
#'
#' # Likelihood profiling
#' \donttest{
#' res <- lik_profile(
#'   x = myscenario,
#'   data = obs,
#'   output = "FrondNo",
#'   par = fit$par,
#'   bounds = list(
#'     k_phot_max = list(0, 30)
#'   ),
#'   refit = FALSE,
#'   type = "fine",
#'   method = "Brent"
#' )
#' # plot
#' plot(res)
#' }
#'
#' @returns A list containing, for each parameter profiled, the likelihood
#' profiling results as a dataframe;
#' the 95% confidence interval; the original parameter value; the likelihood plot object; and
#' the recalibrated parameter values (in case a lower optimum was found)
#' @export

# # License information related to MatLab BYOM v.6.8:
# Copyright (c) 2012-2023 Tjalling Jager
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so

lik_profile <- function(x,
                        par,
                        output,
                        data = NULL,
                        bounds = NULL,
                        refit = TRUE,
                        type = c("coarse", "fine"),
                        individual = FALSE,
                        break_prof = FALSE, # typically we do not want to stop
                        # profiling until we have all parameters, and then make a decision based on parameter space
                        log_scale = FALSE,
                        data_type = c("continuous", "count"),
                        ...) {
  # check inputs
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # convert to list of calibration sets if needed
  if (length(x) == 1) {
    if (is_scenario(x)) {
      if (is.null(data)) {
        stop("Scenario provided, but argument 'data' is missing")
      } else {
        #message("Scenario converted to calibration set")
        if(is.data.frame(data))
          data <- tox_data(data)
        x <- td2cs(x, data, output_var=output)
      }
    }
  }

  # check if type and data_type is one of the 2 options
  type <- match.arg(type)
  data_type <- match.arg(data_type)

  # check correct input parameters
  check_inputs_lik_prof(
    par = par,
    x = x,
    output = output,
    type = type,
    individual = individual,
    data_type = data_type
  )


  # update scenario of calibrationset with the provided par
  for (i in seq_along(x)) {
    x[[i]]@scenario <- x[[i]]@scenario %>%
      set_param(par)
  }

  # initialize parameter list
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # names and boundaries of free parameters
  pnames <- names(par)
  pfree <- list(
    values = par, # free param values
    bounds = get_bounds(x[[1]]@scenario)
  )

  # replace defaults if custom boundaries are given
  if (!is.null(bounds)) {
    check_bounds(bounds)
    pfree$bounds[names(bounds)] <- bounds
  }

  # check that each parameter to profile has boundaries set
  missing <- setdiff(names(par), names(pfree$bounds))
  if(length(missing) > 0) {
    stop("Parameter boundaries missing for parameter", ifelse(length(missing)>1, "s", ""), ": ",
         paste(missing, collapse=", "))
  }

  # Settings for profiling
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # cutoffs for the Chi-square, at a alpha = 0.05 significance level,
  # these values are used in the chi-square test to determine if a fit is
  # significantly different from a previous one (i.e. larger than the cutoff given here)
  chi_crit_j <- stats::qchisq(p=0.95, df=length(pnames))
  chi_crit_s <- stats::qchisq(p=0.95, df=1) # when profiling, we have nested models, differing
  # by 1 parameter (the profiled parameter which is fixed). The likelihood ratio
  # of such nested models follows a Chi-square distribution with 1 degree of freedom


  # setting for profile type
  # values taken from BYOM
  if (type == "fine") {
    # number of iterations
    max_iter <- 50 * (length(pnames) - 1)
    # criteria for change in likelihood ratio
    l_crit_max <- 1 # maximum change in likelihood, above this value, stepsize is decreased
    l_crit_min <- 0.2 # minimum change in likelihood, below this value, stepsize is increased
    l_crit_stop <- chi_crit_s + 5 # stop when likelihood ratio reaches this value
    # settings for step size
    f_step_min <- 0.001 # min stepsize (as fraction of value that is tried)
    f_step_max <- 30 # max stepsize (as fraction of value that is tried)
  } else if (type == "coarse") {
    # number of iterations
    max_iter <- 30 * (length(pnames) - 1)
    # criteria for change in likelihood ratio
    l_crit_max <- 2
    l_crit_min <- 0.7
    l_crit_stop <- chi_crit_s + 2
    # settings for step size
    f_step_min <- 0.005
    f_step_max <- 30
  }


  # The actual profiling
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # a list to store results
  res_list <- list()
  # (important) arguments to this function
  arg_list = list(
    "refit"=refit,
    "type"=type,
    "individual"=individual,
    "log_scale"=log_scale,
    "data_type"=data_type,
    "chi_crit_s"=chi_crit_s
  )
  attr(res_list, "args") <- arg_list

  # profile parameter by parameter
  for (par_select in pnames) {
    # profile the parameter
    res_list[[par_select]] <- suppressWarnings(profile_par(
      par_select = par_select, # par = the parameter to profile
      x = x,
      pfree = pfree,
      type = type,
      output = output,
      max_iter = max_iter,
      f_step_min = f_step_min,
      f_step_max = f_step_max,
      chi_crit_j = chi_crit_j,
      chi_crit_s = chi_crit_s,
      l_crit_max = l_crit_max,
      l_crit_min = l_crit_min,
      l_crit_stop = l_crit_stop,
      refit = refit,
      individual = individual,
      log_scale = log_scale,
      data_type = data_type,
      ...
    ))

    # check if better optimum was found for this parameter (if so, stop profiling)
    if (break_prof) {
      if (any(res_list[[par_select]]$likelihood_profile[, "log_lik_rat"] < 0.00)) {
        # give warning to user
        warning(paste0("Better parameter value found for ", par, ", profiling halted"))
        break
      }
    }
  } # end of parameter for loop

  # Return result list
  class(res_list) <- c(class(res_list), "lik_profile")
  return(res_list)
} # end of lik_profile function



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check inputs for likelihood profiling
#
# @description internal function to check if the parameters of the [lik_profile()]
# function are set correctly
#
# @param x either a single [scenario] or a list of [caliset] objects
# @param par named vector - parameters (names and values) to be profiled
# @param output character vector - the output from the [scenario] or [caliset] that is used in calibration
# @param type "fine" or "coarse" (default) likelihood profiling
# @param individual if 'FALSE' (default), the log likelihood is calculated across
# the whole dataset. Alternatively, if 'TRUE', log likelihoods are calculated for
# each (group of) *set*(s) individually.
# @param data_type Character argument, "continuous" (default) or "count", to specify the data type
# for the log likelihood calculations.
#
# @return x as a list of [caliset], and objects error message when needed

check_inputs_lik_prof <- function(par,
                                  x,
                                  output,
                                  type,
                                  individual,
                                  data_type) {

  # check if attempt to profile more params than possible ~~~~~~~~~~~~~~~~~~~
  if (length(par) > 10) {
    stop("attempt to profile more parameters that function allows, reduce no. of parameters")
  }

  # check model (x)~~~~~~~~~~~~~~~~~~~
  # check if Calibrationset or EffectScenario is provided, and convert EffectScenario to CalibrationSet
  if (is.list(x)) {
    if (is(x[[1]]) != "CalibrationSet") {
      stop("incorrect model specified, please provide a scenario or a list of calibration sets")
    }
  } else { ## input is not a list
    if (is_scenario(x)) {
      # all fine
    } else { # input is not a list, and also not an EffectScenario
      stop("incorrect model specified, please provide a scenario or a list of calibration sets")
    }
  }

  # check if the obs match the data_type
  all <- lapply(x, slot, name = "data")
  observed_data_type <- unlist(lapply(all, function(x) is(x[[2]])))
  stopifnot(any(observed_data_type =="numeric" ))
  if(data_type == "count"){
    stopifnot(any(observed_data_type == "integer"))
  }

  # check par ~~~~~~~~~~~~~~~~~~~
  # check if parameters are provided as named vector
  stopifnot(mode(par) == "numeric")
  if (is.null(names(par))) {
    stop("parameters are not provided as a named vector")
  }

  # check output ~~~~~~~~~~~~~~~~~~~
  stopifnot(mode(output) == "character")


  # check data weights and tags if "individual = TRUE" is chosen ~~~~~~~~~~~~~~~
  if (individual == TRUE) {
    data_weights <- lapply(x, slot, name = "weight")
    data_tag <- unlist(lapply(x, slot, name = "tag"))

    # checks to ensure data weights are given properly
    if (length(data_weights) != length(x)) {
      stop("Please ensure all datasets have a weight assigned to them")
    }

    if (any(unlist(lapply(data_weights, length)) > 1)) {
      stop("Please ensure each individual datasets has only 1 weight assigned to them")
    }

    if (any(is.na(unlist(data_weights)))) {
      stop("Please ensure all datasets have a weight assigned to them (NA is not allowed)")
    }

    # checks to ensure appropriate tags are given
    if (!is.null(data_tag)) {
      for (i in seq_along(unique(data_tag))) {
        index <- which(data_tag == unique(data_tag)[i])
        if (length(unique(unlist(data_weights[index]))) != 1) {
          stop("Please ensure datasets with the same tag, have the same weight assigned to them")
        }
      }
    } else { # if no data_tags are given, make sure that all weights are equal
      if (length(unique(unlist(data_weights))) != 1) {
        stop("Please ensure tag names are provided, when weights are assigned to individual data subsets")
      }
    }
  }
  else {
    data_weights <- as.vector(sapply(x, slot, name="weight"))
    if(length(unique(data_weights)) != 1) {
      warning("Assigned dataset weights will be ignored for profiling if argument `individual=FALSE` ...")
    }
  }
}


# TODO make print statements optional
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter profiling
#
# @description For a specific parameter, this function will sequentially refit
# a model with the parameter value set to
# increasingly lower and higher values, and the likelihood of the model given the
# data calculated (using function [log_lik()]).
#
# The function was inspired by the MatLab BYOM procedure by
# Tjalling Jager, http://debtox.info/byom.html, as described in Jager 2021:
# Robust Likelihood-Based Optimization and Uncertainty Analysis of Toxicokinetic-
# Toxicodynamic Models" Integrated Environmental Assessment and Management
# 17:388-397 DOI: 10.1002/ieam.4333
#
# @param par_select `character`, the parameter to be profiled
# @param x list of [caliset] objects
# @param pfree list of parameter values and their bounds for the profiled parameter
# @param type "fine" or "coarse" (default) likelihood profiling
# @param output `character` vector - the output from the [caliset] that is used during calibration
# @param max_iter `numeric`, maximum number of profiling iterations to attempts
# @param f_step_min, `numeric`,min stepsize (as fraction of value that is tried)
# @param f_step_max, `numeric`,max stepsize (as fraction of value that is tried)
# @param chi_crit_s, `numeric`,cutoffs for the Chi-square under a 95% probability
# @param l_crit_max, `numeric`,maximum change in likelihood, above this value, stepsize is decreased
# @param l_crit_min, `numeric`,minimum change in likelihood, below this value, stepsize is increased
# @param l_crit_stop, `numeric`,stop when likelihood ratio reaches this value
# @param refit if 'TRUE' (default), refit if a better minimum is found
# @param individual if 'FALSE' (default), the log likelihood is calculated across
# the whole dataset. Alternatively, if 'TRUE', log likelihoods are calculated for
# each (group of) *set*(s) individually.
# @param log_scale 'FALSE' (default), option to calculate the log likelihood on a
# log scale (i.e., observations and predictions are log transformed during calculation)
# @param data_type Character argument, "continuous" (default) or "count", to specify the data type
# for the log likelihood calculations.
# @param ... additional parameters passed on to [stats::optim()] and [calibrate()]
#
# @return A list of containing the likelihood profiling results as a dataframe;
# the 95% confidence interval; the original parameter value; the likelihood plot object; and
# the recalibrated parameter values (in case a lower optimum was found)
#' @autoglobal
profile_par <- function(par_select,
                        x, pfree, type, output,
                        max_iter, f_step_min, f_step_max,
                        chi_crit_j, chi_crit_s, l_crit_max, l_crit_min, l_crit_stop,
                        refit, individual, log_scale, data_type, ...) {
  # initialize
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # initialize empty list to store results for this param
  param_res <- list()


  # print name of parameter being profiled
  message("Profiling: ", par_select)

  # get best fit value (i.e., calibrated value)
  x_hat_orig <- pfree$values[[par_select]]

  # get all data, and weights
  all_data <- lapply(x, slot, name = "data")
  data_tag <- unlist(lapply(x, slot, name = "tag"))
  names(all_data) <- data_tag
  data_weights <- unlist(lapply(x, slot, name = "weight"))

  # Make predictions with the original model
  rs <- eval_cs(x, output=output, verbose=FALSE, .ignore_method=TRUE, ...)

  # Calulate log likelihood with original model
  #   Option 1: calculation of loglik across all datasets
  if (individual == FALSE) {
    # calc log lik
    ll_orig <- log_lik(
      npars = length(pfree$values), # Note: only the free params (i.e. ones that you did the calibration on)
      obs = rs$obs,
      pred = rs$pred,
      log_scale = log_scale,
      data_type = data_type
    )
  } else {
    #  Option 2: calculation for individual sub-datasets, which are then combined
    ll_list <- list()
    subsets <- data.frame(obs=rs$obs,
                     pred=rs$pred,
                     wgts=rs$wgts,
                     tags=unlist(rs$tags)) %>%
      dplyr::group_by(tags) %>%
      dplyr::group_split()

    for (i in seq_along(subsets)) {
      ll_list[[i]] <- log_lik(
        npars = length(pfree$values), # Note: only the free params (i.e. ones that you did the calibration on)
        obs = subsets[[i]]$obs,
        pred = subsets[[i]]$pred,
        log_scale = log_scale,
        data_type = data_type
      ) * unique(subsets[[i]]$wgts)
    }
    ll_orig <- sum(unlist(ll_list))
  }


  # additional setting for profile type
  # values taken from BYOM
  if (type == "fine") {
    if (x_hat_orig == 0) { # if parameter value is zero
      f_step_abs <- 1e-4 # This is just a low value that should be ok in most situations
    } else {
      f_step_abs <- 1e-3 * abs(x_hat_orig) # smallest stepsize (in absolute sense)
    }
  } else if (type == "coarse") {
    if (x_hat_orig == 0) {
      f_step_abs <- 1e-3
    } else {
      f_step_abs <- 1e-2 * abs(x_hat_orig)
    }
  }

  message("start param value: ", round(x_hat_orig, 3), "\nLL:", round(ll_orig, 3))

  # refit the model with the par fixed to lower values, and then higher values
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # initialize list to store results during while loop
  log_lik_rat <- c(0)
  par_value <- c(x_hat_orig)
  results <- list(
    loglik = ll_orig,
    log_lik_rat = log_lik_rat,
    par_value = par_value
  )

  # initialize while loop
  reach_limit <- FALSE
  x_try <- NA
  prof_dir <- -1 # first we go down
  change_dir <- FALSE # we are not yet ready to change profiling direction
  iter <- 2 # 2 because first one will be the original value

  while (reach_limit == FALSE) {
    # jump to a lower/upper par value ~~~~~~~~~~~~~~~~~

    if (is.na(x_try)) { # for 1st iteration (and 1st time changing direction), start from the minimum step size (as in BYOM)
      f_step <- f_step_min * abs(x_hat_orig)
      x_try <- x_hat_orig + prof_dir * f_step
    } else { # L. 619 in BYOM
      if (f_step * abs(x_try) >= f_step_abs) { # if stepsize is larger or equal to abs minimum stepsize
        x_try <- x_try + prof_dir * f_step * abs(x_try) # make new x_try value
      } else {
        if (f_step * abs(x_try) < f_step_abs) { # if stepsize is smaller that abs minimum stepsize
          x_try <- x_try + prof_dir * f_step_abs
        } # make new x_try value
        f_step <- f_step_abs / abs(x_try) # define f_step
      }
    }

    # check that values are within the specified parameter bounds
    if (prof_dir == -1) { # while going down
      if (x_try < pfree$bounds[[par_select]][1]) { # if value is below parameter minimum bound
        x_try <- pfree$bounds[[par_select]][[1]] # profile one last time, using the min bound value
        change_dir <- TRUE # then change direction
        warning("Minimum parameter bound reached")
      }
    } else {
      if (prof_dir == 1) { # while going up
        if (x_try > pfree$bounds[[par_select]][2]) { # if value is above parameter maximum bound
          x_try <- pfree$bounds[[par_select]][[2]] # profile one last time, using the max bound value
          reach_limit <- TRUE # stop profiling
          warning("Maximum parameter bound reached")
        }
      }
    }

    # save new param value, x_try
    results[["par_value"]][iter] <- x_try



    # refit model and get likelihood ~~~~~~~~~~~~~~~

    # assign name to x_try
    names(x_try) <- par_select

    # fix the parameter value to the new (lower/higher) value
    for (i in seq_along(x)) {
      x[[i]]@scenario <- x[[i]]@scenario %>% set_param(param = x_try)
    }

    # give initials for remaining free params
    prof_param_index <- which(names(pfree[["values"]]) == par_select)
    pfree_remaining <- unlist(pfree[["values"]])[-prof_param_index]

    # Make predictions with the newly fitted model
    rs <- eval_cs(x, output=output, verbose=FALSE, .ignore_method=TRUE, ...)

    # Calulate log likelihood with the newly fitted model
    #   Option 1: calculation of loglik across all datasets
    if (individual == FALSE) {
      # calc log lik
      ll_new <- log_lik( # list with only 1 entry
        npars = length(pfree$values), # Note: only the free params (i.e. ones that you did the calibration on)
        obs = rs$obs,
        pred = rs$pred,
        log_scale = log_scale,
        data_type = data_type
      )
    } else {
      #  Option 2: calculation for individual sub-datasets, which are then combined
      ll_list <- list()
      subsets <- data.frame(obs=rs$obs,
                            pred=rs$pred,
                            wgts=rs$wgts,
                            tags=unlist(rs$tags)) %>%
        dplyr::group_by(tags) %>%
        dplyr::group_split()

      for (i in seq_along(subsets)) {
        ll_list[[i]] <- log_lik(
          npars = length(pfree$values), # Note: only the free params (i.e. ones that you did the calibration on)
          obs = subsets[[i]]$obs,
          pred = subsets[[i]]$pred,
          log_scale = log_scale,
          data_type = data_type
        ) * unique(subsets[[i]]$wgts)
      }
      ll_new <- sum(unlist(ll_list))
    }

    # save new LL
    results[["loglik"]][iter] <- ll_new

    # save -2*ln(likelihood ratio) comparing the full (original) model with the nested (1 param fixed) model
    results[["log_lik_rat"]][iter] <- -2 * (ll_new - ll_orig) # matlab BYOM L. 704

    # difference with ll from previous iteration (needed to decide on step size etc.) # matlab BYOM L. 705
    if (iter == 1) {
      delta_l <- abs(results[["loglik"]][iter] - ll_orig)
    } else {
      delta_l <- abs(results[["loglik"]][iter] - results[["loglik"]][iter - 1])
    }

    # evaluate likelihood ~~~~~~~~~~~~~~~~~~~

    # X2 difference test to compare nested models
    if (results[["log_lik_rat"]][iter] > chi_crit_s) { # if the point is outside the 95% conf. region
      if (prof_dir == -1) { # while going down
        change_dir <- TRUE # then change direction
        message("hit 95% Conf Region, changing direction")
      } else {
        if (prof_dir == 1) { # while going up
          reach_limit <- TRUE # stop profiling
          message("hit 95% Conf Region, end of profiling")
        }
      }
    }

    if (delta_l > l_crit_stop) { # stop while-loop if probability is high enough
      reach_limit <- TRUE
      message("critical limit reached, end of profiling")
    } else {
      if (delta_l < l_crit_min) { # if change in likelihood is very small
        f_step <- min(2 * f_step, f_step_max) # double the stepsize, unless greater than max step
      }
    }


    # if needed, change direction
    if (change_dir == TRUE) {
      prof_dir <- 1 # change direction
      x_try <- x_hat_orig # restart from original param value
    }

    iter <- iter + 1
  } # end of while loop


  # reform results list into dataframe
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  res_df <- data.frame(
    par_value = results[["par_value"]],
    log_lik_rat = results[["log_lik_rat"]],
    loglik = results[["loglik"]]
  )


  # get prof region
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  min <- min(res_df$par_value)
  max <- max(res_df$par_value)

  # get 95% CI
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # indicate the first value
  res_df$start <- "No"
  res_df[1, "start"] <- "Yes"
  res_df$start <- as.factor(res_df$start)

  # order dataset, from low to high param value
  res_df <- res_df[order(res_df$par_value), ]

  orig_ind <- which(res_df$start == "Yes")
  # lower 95
  low_df <- res_df[1:orig_ind, ]
  if (min(low_df$log_lik_rat) < 0) { # in case a lower minimum is found within the data subset
    min_ind <- which(low_df$log_lik_rat == min(low_df$log_lik_rat))
    low_df <- low_df[1:min_ind, ]
  }
  if (nrow(low_df) > 1 && !all(low_df[, "log_lik_rat"] == 0)) {
    # needs to be at least 2 values, at least one with a ll not 0
    f_inter <- stats::approxfun(
      y = low_df$par_value,
      x = low_df$log_lik_rat, rule = 2
    )
    low95 <- f_inter(chi_crit_s)
  } else {
    low95 <- low_df$par_value
    if (length(low95) > 1) {
      low95 <- min(low95)
    } # in case multiple values with ll of 0, take the min param value
  }
  # upper 95
  up_df <- res_df[orig_ind:nrow(res_df), ]
  if (min(up_df$log_lik_rat) < 0) { # in case a lower minimum is found within the data subset
    min_ind <- which(up_df$log_lik_rat == min(up_df$log_lik_rat))
    up_df <- up_df[min_ind:nrow(up_df), ]
  }
  if (nrow(up_df) > 1 && !all(up_df[, "log_lik_rat"] == 0)) {
    # needs to be at least 2 values, at least one with a ll not 0
    f_inter <- stats::approxfun(
      y = up_df$par_value,
      x = up_df$log_lik_rat, rule = 2
    )
    up95 <- f_inter(chi_crit_s)
  } else {
    up95 <- up_df$par_value
    if (length(up95) > 1) {
      up95 <- max(up95)
    } # in case multiple values with ll of 0, take the max param value
  }

  # save all outputs for the parameter
  param_res <- list(
    likelihood_profile = res_df,
    confidence_interval = c(low95, up95),
    prof_region = c(min, max),
    orig_par_value = x_hat_orig,
    fit_new = NULL
  )

  # recalibrate if needed
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Check in results if any log likelihood ratio is lower than the original
  if (refit == TRUE) {
    if (any(res_df[["loglik"]] > ll_orig + 0.01)) {
      # only if the difference is at least 0.01 as in BYOM (value taken from BYOM)
      # this ignores tiny, irrelevant, improvement
      best_par_value <- res_df[which(res_df$log_lik_rat == min(res_df$log_lik_rat)), "par_value"]
      names(best_par_value) <- par_select
      # update par value in model
      for (i in seq_along(x)) {
        x[[i]]@scenario <- x[[i]]@scenario %>% set_param(param = best_par_value)
      }
      # recalibrate
      fit_new <- calibrate(
        x = x,
        par = pfree_remaining,
        output = output,
        log_scale = log_scale,
        verbose = FALSE,
        ...
      )
      # save recalibrated result
      param_res[["fit_new"]] <- list(
        best_fit_param = c(best_par_value, fit_new$par),
        recalibration_outputs = fit_new
      )
    } else { # if no better optimum was found for this parameter
      param_res[["fit_new"]] <- "No recalibration was necessary."
    }
  }
  return(param_res)
} # end of profile_par function


#' Calculate log likelihood
#'
#' @description Calculates the sum of log likelihoods of each observation given
#' the model parameterization.
#'
#' Current implementations enable log likelihood calculations for:
#' 1) continuous data, considering a normal distribution around the prediction
#' for each datapoint,
#' 2) count data, considering a multinomial distribution for data reporting the
#' number of survivors over time.
#'
#' The log likelihood calculation for count data was inspired by a MatLab BYOM v.6.8
#' procedure, created by Tjalling Jager. For details, please refer to BYOM
#' (http://debtox.info/byom.html) as well as Jager (2021).
#'
#' @references
#' Jager T, 2021. Robust Likelihood-Based Optimization and Uncertainty Analysis
#' of Toxicokinetic-Toxicodynamic Models. Integrated Environmental Assessment and
#' Management 17:388-397. \doi{10.1002/ieam.4333}
#'
#' @param obs numeric vector of observed values
#' @param pred numeric vector of predicted values
#' @param data_type determines the if likelihood profiling is conducted for
#' `"continuous"` (default) or `"count"` data
#' @param log_scale `FALSE` (default), option to calculate the log likelihood on a
#' log scale (i.e., observations and predictions are log transformed during calculation)
#' @param npars named numeric vector of parameters that the model was calibrated on,
#' required for `"continuous"` data type, optional for `"count"`.
#'
#' @return the log likelihood value
#' @export
#'
#' @examples
#' # simple example for continuous data #####
#' # observations
#' obs <- c(12, 38, 92, 176, 176, 627, 1283, 2640)
#' # intercept, a, and slope, b, of a Poisson regression fitted through obs
#' pars <- c(a = 2, b = 0.73)
#' # predictions with the Poisson regression
#' pred <- c(15.43, 32.15, 66.99, 139.57, 290.82, 605.94, 1262.52, 2630.58)
#' # example plot
#' plot(seq(1:length(obs)), obs)
#' lines(seq(1:length(obs)), pred)
#' log_lik(
#'   obs = obs,
#'   pred = pred,
#'   npars = length(pars),
#' )
#'
#' # example with count data and GUTS model #####
#' library(dplyr)
#' # observational data
#' dt <- ringtest_c %>% filter(replicate == "E")
#' myexposure <- dt %>% select(time, conc)
#' obs <- dt %>%
#'   mutate(S=Nsurv / max(Nsurv)) %>%
#'   select(time, S)
#' # GUTS model
#' GUTS_RED_IT() %>%
#'   set_param(c(hb = 0)) %>%
#'   set_exposure(myexposure) -> myscenario
#' # fit
#' fit <- calibrate(
#'   x = myscenario,
#'   par = c(kd=1.2, alpha=9.2, beta=4.3),
#'   data = obs,
#'   output = "S")
#' # update
#' myscenario <- myscenario %>% set_param(fit$par)
#' # simulate
#' pred <- myscenario %>% simulate()
#' pred <- pred$S #* max(obs$S)
#' obs <- obs$S
#' # calc likelihood
#' log_lik(obs,
#'     pred,
#'     data_type = "count")
#'

log_lik <- function(obs,
                    pred,
                    data_type = c("continuous", "count"),
                    log_scale = FALSE,
                    npars = NULL){

  # general checks
  stopifnot(length(obs) == length(pred))
  if(any(is.na(c(obs, pred)))){
    warning("missing values in observations or predictions")
  }
  data_type <- match.arg(data_type)

  # transformations
  if(log_scale == TRUE){
    obs = log(obs)
    pred = log(pred)
  }

  # Likelihoods
  if(data_type == "continuous"){

    # specific checks
    stopifnot(is.numeric(npars))
    stopifnot(length(npars) == 1)

    # calculations
    k <- npars
    res <- obs - pred
    n <- length(res)
    SSE <- sum(res^2)
    sigma <- sqrt(SSE / (n - k))
    sigma_unbiased <- sigma * sqrt((n - k) / n)

    # log likelihood for normal probability density
    LL <- sum(log(stats::dnorm(x = obs, mean = pred, sd = sigma_unbiased)))
    return(LL)
  }
  if(data_type == "count"){

    # calculations
    Ndeaths <- -diff(obs)       # numbers of deaths
    Mdeaths <- -diff(pred)      # conditional probabilities of deaths
    Mdeaths <- max(Mdeaths,1e-50)   # otherwise we get problems taking the logarithm
    pred    <- max(pred,1e-50)      # otherwise we get problems taking the logarithm

    # log likelihood for count data (multinomial)
    LL  <- sum( Ndeaths * log(Mdeaths) )
    return(LL)
  }

}

