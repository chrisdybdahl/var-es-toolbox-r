require(rugarch)
require(xts)
require(zoo)

forecast_u_GARCH <- function(
  data,
  c,
  n,
  m,
  p = 1,
  q = 1,
  r = 1,
  model = "sGARCH",
  dist = "norm",
  solver = "hybrid",
  solver.control =list(tol = 1e-7),
  ...
) {
  df <- tail(data, n + m)
  data_xts <- xts::xts(df$Return, order.by = df$Date)

  spec <- rugarch::ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(garchOrder = c(p, q), model = model),
    distribution.model = dist
  )

  garch_roll <- rugarch::ugarchroll(
    spec = spec,
    data = data_xts,
    forecast.length = n,
    window.size = m,
    refit.every = r,
    refit.window = 'moving',
    solver = solver,
    solver.control = solver.control,
    ...
  )

  sigmaFor <- garch_roll@forecast$density$Sigma

  # Assume mean is zero
  # TODO: Assume mean is zero
  if (dist == "norm") {
    quantile <- sapply(c, function(conf) stats::qnorm(conf))
    pdf_value <- sapply(c, function(conf) stats::dnorm(quantile[which(c == conf)]))
    quantile <- matrix(rep(quantile, n), nrow = n, byrow = TRUE)
    pdf_value <- matrix(rep(pdf_value, n), nrow = n, byrow = TRUE)
    var <- -sigmaFor * (quantile)
    es <- sigmaFor * (pdf_value / c)
  } else if (dist == "std") {
    shapeFor <- garch_roll@forecast$density$Shape
    quantile <- sapply(c, function(conf) rugarch::qdist("std", p = conf, mu = 0, sigma = 1, shape = shapeFor))
    pdf_value <- sapply(c, function(conf) rugarch::ddist("std", quantile[which(c == conf)], mu = 0, sigma = 1, shape = shapeFor))
    var <- -sigmaFor * (quantile)
    es <- sigmaFor * ((shapeFor + quantile^2) / (shapeFor - 1) * (pdf_value / c))
  } else if (dist == "sstd") {
    shapeFor <- garch_roll@forecast$density$Shape
    skewFor <- garch_roll@forecast$density$Skew
    quantile <- sapply(c, function(conf) rugarch::qdist("sstd", p = conf, mu = 0, sigma = 1, skew = skewFor, shape = shapeFor))
    pdf_value <- sapply(c, function(conf) rugarch::ddist("sstd", quantile[which(c == conf)], mu = 0, sigma = 1, skew = skewFor, shape = shapeFor))
    var <- -sigmaFor * (quantile)
    es <- sigmaFor * ((shapeFor + quantile^2) / (shapeFor - 1) * (pdf_value / c))
  } else {
    stop("Unsupported distribution type. Use 'norm', 'std', or 'sstd'.")
  }

  dates <- tail(data$Date, n)
  VaR <- xts::xts(var, order.by = dates)
  ES <- xts::xts(es, order.by = dates)
  VOL <- xts::xts(as.numeric(sigmaFor), order.by = dates)
  colnames(VaR) <- paste0("VaR_", c)
  colnames(ES) <- paste0("ES_", c)
  colnames(VOL) <- "VOL"
  results_xts <- xts::merge.xts(VaR, ES, VOL)

  if (any(is.na(results_xts))) {
    warning("There were NA values, carry forward last non-NA")
    results_xts <- zoo::na.locf(results_xts, na.rm = FALSE)
  }

  return(results_xts)
}
