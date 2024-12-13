require(xts)
require(caviar)

forecast_u_CAViaR <- function(data, c, n, m, verbose = 3, ...) {
  # Initialize empty matrices to store results
  var_matrix <- matrix(NA, nrow = n, ncol = length(c))
  es_matrix <- matrix(NA, nrow = n, ncol = length(c))

  # Iterate over each confidence level
  for (i in seq_along(c)) {
    # Call the original forecast_u_CAViaR function for each confidence level
    result <- caviar::RollCAViaR(data, c[i], n, m, verbose = verbose, control.optim = list(factr = 1e-7), control.solnp = list(tol = 1e-7, trace = 0), ...)
    print(result)
    # Extract VaR and ES for the current confidence level
    var_matrix[, i] <- as.numeric(result$VaR)
    es_matrix[, i] <- as.numeric(result$ES)

    print(tail(var_matrix, 10))
  }

  dates <- tail(data$Date, n)
  VaR_xts <- xts::xts(var_matrix, order.by = dates)
  ES_xts <- xts::xts(es_matrix, order.by = dates)
  colnames(VaR_xts) <- paste0("VaR_", c)
  colnames(ES_xts) <- paste0("ES_", c)
  results_xts <- merge(VaR_xts, ES_xts)

  print("done")

  return(results_xts)
}
