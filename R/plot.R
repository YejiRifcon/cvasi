#' S3 plotting functions
#'
#' These functions overload [base::plot()] to provide simple plotting
#' routines to display various time-series and scenario objects.
#'
#' @name plot
#' @param x object to plot
#' @param y unused parameter
#' @param ... unused parameters
#' @returns (*ggplot2*) plot object
NULL

# S3 overload to plot dose response curves
#
#' @param scale_x character, controls how the x-axis is scaled. `log10` for
#'   a log10-scaled axis, `none` for no scaling, and `auto` for automatic selection
#' @describeIn plot Plot dose response curves
#' @export
plot.cvasi_drc <- function(x, y, scale_x=c("auto", "log10", "none"), ...) {
  if(!missing(y))
    warning("Parameter `y` is unused, ignoring argument")
  if(!missing(...))
    warning("Parameter `...` is unused, ignoring additional arguments")

  scale_x <- match.arg(scale_x)
  if(scale_x == "auto") {
    exp_min <- floor(log10(min(x$mf)))
    exp_max <- ceiling(log10(max(x$mf)))
    if(is.na(exp_min) | is.na(exp_max))
      scale_x <- "none"
    else if(abs(exp_max - exp_min) > 2)
      scale_x <- "log10"
    else
      scale_x <- "none"
  }

  plot <- ggplot2::ggplot(x) +
    ggplot2::geom_line(ggplot2::aes(mf, effect, color=endpoint))
  if(scale_x == "log10") {
    plot <- plot + ggplot2::scale_x_log10()
  }

  # create a nice plot title
  subtitle <- ggplot2::element_blank()
  sc <- attr(x, "scenario")
  if(!is.null(sc)) {
    nm <- get_model_name(sc)
    tg <- get_tag(sc)
    if(is.na(tg))
      subtitle <- nm
    else
      subtitle <- paste0(nm, " #", tg)
  }

  plot +
    ggplot2::scale_color_discrete(name="Endpoint") +
    ggplot2::labs(title="Dose response curve",
                  subtitle=subtitle,
                  x=ifelse(scale_x == "log10", "Multiplication factor [log-scale] (-)", "Multiplication factor (-)"),
                  y="Effect (-)") +
    ggplot2::theme_bw()
}

#' @describeIn plot Plot return value of [simulate()]
#' @export
plot.cvasi_simulate <- function(x, y, ...) {
  if(!missing(y))
    warning("Parameter `y` is unsused, ignoring argument")
  if(!missing(...))
    warning("Parameter `...` is unused, ignoring additional arguments")

  # pivot table in wide-format to long-format. we have to assume that the
  # first column represents time. the remaining numerical columns will be
  # plotted.
  tidyr::pivot_longer(x, cols=!c(1) & dplyr::where(is.numeric)) %>%
    dplyr::rename(time=1) %>%                                  # make sure we have a well defined name
    dplyr::mutate(name=factor(name, levels=unique(name))) %>%  # enforce ordering of state variables in plot
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(time, value, color=name)) +
    ggplot2::facet_wrap(~ name, scales="free") +
    ggplot2::labs(title="Simulation result", x="Time", y="Value (?)") +
    ggplot2::guides(color="none") +
    ggplot2::theme_bw()
}

#' @describeIn plot Plot likelihood profiles.
#' @autoglobal
#' @export
plot.lik_profile <- function(x, y, ...) {
  if(!missing(y))
    warning("Parameter `y` is unsused, ignoring argument")
  if(!missing(...))
    warning("Parameter `...` is unused, ignoring additional arguments")

  npars <- length(x)
  n1 <- ceiling(sqrt(npars))

  args <- attr(x, "args", exact=TRUE)

  # create one plot per profiled parameter
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  plots <- list()
  for(nm in names(x)) {
    df <- x[[nm]][["likelihood_profile"]]
    plot <- df %>%
      ggplot2::ggplot() +
      # # zone of acceptabiity
      # ggplot2::geom_rect(
      #   ggplot2::aes(xmin=-Inf, xmax=Inf, ymin=min(log_lik_rat), ymax=args$chi_crit_s),
      #   fill = "darkolivegreen3", alpha = 0.2
      # ) +
      ggplot2::annotate("rect", xmin=-Inf, xmax=Inf, ymin=min(df$log_lik_rat), ymax=args$chi_crit_s,
                        alpha=0.3, fill="darkolivegreen3") +
      ggplot2::geom_hline(yintercept=args$chi_crit_s, color="tomato", linewidth=1.2) +
      # points for values
      ggplot2::geom_point(ggplot2::aes(x=par_value, y=log_lik_rat, color=start), size=2) +
      ggplot2::scale_color_manual(values = c("black", "orange")) +
      # connect points
      ggplot2::geom_line(ggplot2::aes(x=par_value, y=log_lik_rat)) +
      # could add 95% CI
      ggplot2::geom_vline(xintercept=x[[nm]][["confidence_interval"]], linetype="dashed", color="grey30", linewidth=.5) +
      ggplot2::labs(x=nm, y="LL ratio (-)") +
      # cosmetics
      ggplot2::scale_y_continuous(limits=c(min(df$log_lik_rat), NA), expand=c(0, 0)) +
      ggplot2::theme_classic()
    plots <- append(plots, list(plot))
  }
  do.call(gridExtra::grid.arrange, c(plots, list(nrow = n1)))
}

#' @describeIn plot Alias of `plot.lik_profile()` for backwards-compatibility.
#' @export
plot_lik_profile <- function(x) {
  lifecycle::deprecate_soft("1.5.0", "cvasi::plot_lik_profile()", "plot()")
  plot(x=x)
}


#' Plot profiled parameter space
#'
#' The function provides bivariate parameter space plots indicating
#' parameter draws (from the 95% confidence intervals per parameter obtained
#' through likelihood profiling) that fall within the inner rim (in green, i.e.
#' parameter sets which are not significantly different from the original, based
#' on a chi-square test). The original parameter set is also indicated (in orange),
#' and, if different from the original set, the best fit parameter set is indicated
#' (in red).
#'
#' @param x Return value of [explore_space()]
#' @param y Unused argument
#' @param ... Unused arguments
#' @returns plot object
#' @autoglobal
#' @export
plot.param_space <- function(x, y, ...){
  if(!missing(y))
    warning("Parameter `y` is unsused, ignoring argument")
  if(!missing(...))
    warning("Parameter `...` is unused, ignoring additional arguments")

  # find best LLR value, and compare with original LLR (which is per definition = 0)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (min(x$LLR) < 0) {
    plot_colors <- c("tomato", "orange", "#4DAC26",  "#0571B0")
    plot_size <- c(3, 3, 1.5, 1.5)
    index <- which(x$LLR_quality %in% c("Original fit", "Best fit"))
    x$LLR_quality <- as.factor(x$LLR_quality)
    x$LLR_quality <- factor(x$LLR_quality,
                            levels = c("Best fit", "Original fit", "Inner", "Outer"))
  } else {
    plot_colors <- c( "orange","#4DAC26", "#0571B0")
    plot_size <- c(3, 1.5, 1.5)
    index <- which(x$LLR_quality %in% "Original fit")
    x$LLR_quality <- as.factor(x$LLR_quality)
    x$LLR_quality <- factor(x$LLR_quality,
                            levels = c("Original fit", "Inner", "Outer"))
  }


  # reorder so that the original and best fit are the last datapoints (and
  # thus are plotted on top)
  x_best <- x[index,]
  x_withoutbest <- x[-index,]
  x <- rbind(x_withoutbest, x_best)

  # plot of all correlations
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # dummy plot to get a legend
  p <- x %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = colnames(x)[1],
        y = LL,
        color = LLR_quality)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = plot_colors)
  leg <- GGally::grab_legend(p)

  # actual plot, complete
  full_plot <- x %>%
    GGally::ggpairs(
      columns = colnames(x)[-which(colnames(x) %in% c("LL", "LLR", "LLR_quality"))],
      ggplot2::aes(
        color = LLR_quality,
        size = as.factor(LLR_quality)
      ),
      upper = list(continuous = "blank"),
      diag = list(continuous = "blankDiag"),
      title = "Parameter space",
      legend = leg
    ) +
    ggplot2::scale_color_manual(values = plot_colors) +
    ggplot2::scale_size_manual(values = plot_size) +
    ggplot2::theme_classic()

  # function to modify ggpairs output (trim empty space)
  gpairs_lower <- function(g) {
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow - 1

    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1

    g
  }

  # modification of the plot
  full_plot <- gpairs_lower(full_plot)

  return(full_plot)
}

#' @describeIn plot.param_space Alias of `plot.param_space()` for backwards-compatibility.
#' @export
plot_param_space <- function(x) {
  lifecycle::deprecate_soft("1.5.0", "cvasi::plot_param_space()", "plot()")
  plot(x=x)
}
