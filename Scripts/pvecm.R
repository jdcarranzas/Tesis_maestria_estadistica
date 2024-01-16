require(pvecm)

pvecm(
  data[c(6,7,8)],
  data[c(9)],
  I_0_dols = NULL,
  cross_sections = data[c(5)],
  time = data[c(3)],
  dummies = NULL,
  method = "FM",
  deterministic_long = "none",
  deterministic_short = "none",
  vecm_lags = 1,
  maximum_lags = 3,
  kernel = "ba",
  aic_small = TRUE,
  bandwidth = "and",
  n.lead = NULL,
  n.lag = NULL,
  kmax = "k4",
  info.crit = "AIC",
  demeaning = FALSE,
  check = TRUE
)
