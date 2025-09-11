
# data
ringtest_c <- read.table(file = "data-raw/efsa_ringtest_c.txt",
                              header = TRUE) %>%
  dplyr::select(time, Nsurv, replicate, conc)

# write data
usethis::use_data(ringtest_c, overwrite=TRUE)

rm(ringtest_c)
