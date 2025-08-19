# Example scenario using the algae model by Weber et al. (2012) for R. subcapitata
# exposed to isoproturon in reactor A. Model parameters are set according to the
# values reported by Weber et al. Exact values of the exposure series were not
# available, but it was reconstructed from information given in text and figures.
#
# Parameter values reported by Weber et al. appear to be rounded, so small
# deviations in simulation results may be possible when compared to published
# figures etc.
#
# Please note: the algae species is named *R. subcapitata* in EFSA (2018),
#   and *P. subcapitata* (older name) in Weber et al.
#
# @references
# Weber D, Schaefer D, Dorgerloh M, Bruns E, Goerlitz G, Hammel K, Preuss TG
# and Ratte HT, 2012. Combination of a higher-tier flow-through system and
# population modeling to assess the effects of time-variable exposure of
# isoproturon on the green algae Desmodesmus subspictatus and
# Pseudokirchneriella subcapitata. Environmental Toxicology and
# Chemistry, 31(4), 899-908. \doi{10.1002/etc.1765}


# Reconstructed exposure time-series for time-frame from day zero to 42,
# cf. Figures 2A and 4A in Weber et al. (2012)
times <- seq(0, 42, 0.1)
conc <- c()
D <- 0.5 # dilution rate in reactor (1/day)

# The dilution of substance in the reactor by flow-through is modeled by
# using exponential decay with a half-live equal to `D`, i.e. DT50=0.5 days
for(t in times) {
  c <- 0
  if(t >= 13 & t < 15)
    c <- 65 * exp(D * (-t + 13))
  else if(t >= 15 & t < 23)
    c <- 160 * exp(D * (-t + 15))
  else if(t >= 23 & t < 30)
    c <- 160
  else if(t >= 30)
    c <- 160 * exp(D * (-t + 30))
  conc <- c(conc, c)
}
weber_exposure <- data.frame(t=times, c=conc)
rm(times, conc, D, c, t)

# Initial state,
# cf. Table 2 in Weber et al. (2012)
y0 <- c(A = 100, Q = 0.01, P = 0.36 * 0.5)

# Model parameters for R. subcapitata with median EC50,
# cf. Tables 2 and C in Weber et al. (2012)
params <- c(mu_max = 1.7820, m_max = 0.0500, v_max = 0.0620, k_s = 0.0625,
            Q_min = 0.0011, Q_max = 0.0175,
            T_opt = 28, T_min = 0, T_max = 42, I_opt = 120,
            EC_50 = 128, b = 1.199
)

# Run a control simulation to find steady-state of algae biomass
control <- Algae_Weber() %>%
  set_init(y0) %>%
  set_param(params) %>%
  set_noexposure() %>%
  set_forcings(I = 100, T_act = 24) %>%
  set_times(seq(0, 72))
y1 <- control %>%
  simulate() %>%
  tail(n=1) %>%
  dplyr::select(A, Q, P) %>%
  unlist()

# Create scenario that starts in steady-state and reproduces conditions
# of R. subcapitata exposed to isoproturon in reactor A, cf. Figure 4A in
# Weber et al. (2012)
rsubcapitata <- control %>%
  set_init(y1) %>%
  set_exposure(weber_exposure) %>%
  set_times(seq(0, 42, 0.1))

usethis::use_data(rsubcapitata, overwrite=TRUE)

rm(weber_exposure, control, rsubcapitata, params, y0, y1)




