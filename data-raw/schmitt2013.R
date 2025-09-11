# Create rda file for Schmitt (2013) lemna effect data containing multiple
# treatments

#' Schmitt W., Bruns E., Dollinger M., and Sowig P., 2013: *Mechanistic TK/TD-model
#' simulating the effect of growth inhibitors on Lemna populations*. Ecol Model 255,
#' pp. 1-10. \doi{10.1016/j.ecolmodel.2013.01.017}

# ------------------------------------------------------------------------------
# Exposure data
schmitt2013 <- read.table(file="data-raw/schmitt2013.txt", header=TRUE, sep="\t") %>%
  dplyr::select(time=t, obs, trial=ID, conc) %>%
  # values at 7.01 aren't real observations, but just a timepoint inserted for
  # technical purposes, i.e. to represent the exposure step-function
  dplyr::mutate(obs=ifelse(time == 7.01, NA_real_, obs))

usethis::use_data(schmitt2013, overwrite=TRUE)

rm(schmitt2013)
