
# data
ringtest_sd <- read.table(file = "data-raw/efsa_ringtest_a_SD.txt",
                              header = TRUE) %>%
  dplyr::select(time, Nsurv, replicate, conc)

# write data
usethis::use_data(ringtest_sd, overwrite=TRUE)

rm(ringtest_sd)


# data
ringtest_it <- read.table(file = "data-raw/efsa_ringtest_a_IT.txt",
                                 header = TRUE) %>%
  dplyr::select(time, Nsurv, replicate, conc)

# write data
usethis::use_data(ringtest_it, overwrite=TRUE)

rm(ringtest_it)
