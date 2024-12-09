require(GAS)
require(rugarch)
require(xts)
require(zoo)

forecast_u_GAS <- function(
    data,
    c,
    n,
    m,
    r = 1,
    dist = "norm",
    ...
) {
  df <- tail(data$Return, n + m)
  GASSpec <- GAS::UniGASSpec(
    Dist = dist,
    ScalingType = "Identity",
    GASPar = list(scale = TRUE)
  )

  gas_roll <- GAS::UniGASRoll(
    data = df,
    GASSpec = GASSpec,
    ForecastLength = n,
    RefitEvery = r,
    RefitWindow = "moving",
    ...
  )

  var <- -GAS::quantile(gas_roll, probs = c)
  es <- -GAS::ES(gas_roll, probs = c)

  # Create xts objects for VaR and ES
  dates <- tail(data$Date, n)
  VaR <- xts::xts(var, order.by = dates)
  ES <- xts::xts(es, order.by = dates)
  colnames(VaR) <- paste0("VaR_", c)
  colnames(ES) <- paste0("ES_", c)
  results_xts <- merge(VaR, ES)

  if (any(is.na(results_xts))) {
    warning("There were NA values, carry forward last non-NA")
    results_xts <- zoo::na.locf(results_xts, na.rm = FALSE)
  }

  return(results_xts)
}
