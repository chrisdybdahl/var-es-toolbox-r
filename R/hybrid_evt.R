require(rugarch)
require(xts)
require(extRemes)
require(zoo)
require(progress)

forecast_u_EVT_GARCH <- function(
  data,
  c,
  n,
  m,
  p = 1,
  q = 1,
  r = 1,
  t = 0.9,
  verbose = 1,
  model = "sGARCH",
  dist = "norm",
  min_exceedances = 10,
  solver = "hybrid",
  solver.control = list(trace = 0, tol = 1e-7),
  ...
) {

  if (verbose > 0) {
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsed || Estimated time remaining: :eta]",
                                     total = n,
                                     complete = "=",
                                     incomplete = "-",
                                     current = ">",
                                     show_after = 0.2,
                                     clear = TRUE)
  }

  df <- tail(data, n + m)
  data_xts <- xts::xts(df$Return, order.by = df$Date)

  spec_unrestricted <- rugarch::ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    variance.model = list(garchOrder = c(p, q), model = model),
    distribution.model = dist
  )

   # Initialize storage for results
  var <- matrix(NA, nrow = n, ncol = length(c), dimnames = list(NULL, paste0("VaR_", c)))
  es <- matrix(NA, nrow = n, ncol = length(c), dimnames = list(NULL, paste0("ES_", c)))

  vol <- numeric(n)
  xi <- NA
  beta <- NA
  params <- NA
  for (i in 1:n) {
    if (verbose > 0) pb$tick()
    window_start <- i
    window_end <- m + i - 1
    window <- data_xts[window_start:window_end]

    evt_adj <- FALSE
    valid_out <- FALSE

    # Reparametrize for every r period
    if (i %% r == 0 || i == 1) {
      fit <- rugarch::ugarchfit(
        spec = spec_unrestricted,
        data = df,
        solver = solver,
        solver.control = solver.control
      )

      params <- rugarch::coef(fit)

      std_res <- rugarch::residuals(fit, standardize = TRUE)

      # Threshold selection
      u <- quantile(std_res, probs = t)

      # Calculate exceedances with a marked point process framework
      n_u <- length(std_res[std_res > u])

      # Fit the Generalized Pareto Distribution with tryCatch
      evt_fit <- tryCatch({
        extRemes::fevd(std_res, threshold = u, type = "GP", method = "MLE")
      }, error = function(e) {
        warning("GPD fitting failed: ", e$message)
        return(NULL)
      })

      if (!is.null(evt_fit)) {
        # Extract parameters if fitting was successful
        xi <- evt_fit$results$par["shape"]
        beta <- evt_fit$results$par["scale"]
      } else {
        # Handle failed fitting
        warning("GPD fitting was not successful on date: ", as.character(df$Date[m + i - 1]), ". Skipping EVT adjustment.")
      }
    }

    # day ahead
    spec_restricted <- rugarch::ugarchspec(
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      variance.model = list(garchOrder = c(p, q), model = model),
      distribution.model = dist,
      fixed.pars = list(params)
    )

    fit <- rugarch::ugarchfit(
      spec = spec_restricted,
      #spec = spec_unrestricted,
      data = window,
      solver = solver,
      solver.control = solver.control,
      ...
    )

    vol[i] <- rugarch::sigma(ugarchforecast(fit, n.ahead = 1))
    #

    # current day
    #vol[i] <- rugarch::sigma(fit)[m]
    #


    if (!any(is.na(xi)) && !any(is.na(beta)) && (n_u != 0)) {
      # Using Pickands-Balkema-de Haan Extreme Value Theorem (Balkema & de Haan, 1974) we can derive VaR for GPD
      # From Hull (2018) we have a definition for ES for GPD

      var_evt <- if (xi != 0) {
        abs(u) + (beta / xi) * ((1 - c)^(-xi) - 1)
      } else {
        abs(u) - beta * log(1 - c)
      }

      es_evt <- if (xi != 0) {
        (var_evt + beta - abs(u) * xi) / (1 - xi)
      } else {
        var_evt + beta
      }

      # Ensure the calculated values are valid
      if (any(is.na(var_evt)) || any(is.na(es_evt)) || any(var_evt <= 0) || any(es_evt <= 0)) {
        warning("VaR or ES calculation returned NA or non-positive values on date: ", as.character(df$Date[m + i - 1]), ", fallback to GARCH")

        sigmaFor <- vol[i]

        if (dist == "norm") {
          quantile <- sapply(c, function(conf) stats::qnorm(conf))
          pdf_value <- sapply(c, function(conf) stats::dnorm(quantile[which(c == conf)]))
          var[i, ] <- -sigmaFor * (quantile)
          es[i, ] <- sigmaFor * (pdf_value / c)
        } else if (dist == "std") {
          shapeFor <- fit@fit$coef["shape"]
          quantile <- sapply(c, function(conf) rugarch::qdist("std", p = conf, mu = 0, sigma = 1, shape = shapeFor))
          pdf_value <- sapply(c, function(conf) rugarch::ddist("std", quantile[which(c == conf)], mu = 0, sigma = 1, shape = shapeFor))
          var[i, ] <- -sigmaFor * (quantile)
          print(var[i, ])
          es[i, ] <- sigmaFor * ((shapeFor + quantile^2) / (shapeFor - 1) * (pdf_value / c))
        } else if (dist == "sstd") {
          shapeFor <- fit@fit$coef["shape"]
          skewFor <- fit@fit$coef["skew"]
          quantile <- sapply(c, function(conf) rugarch::qdist("sstd", p = conf, mu = 0, sigma = 1, skew = skewFor, shape = shapeFor))
          pdf_value <- sapply(c, function(conf) rugarch::ddist("sstd", quantile[which(c == conf)], mu = 0, sigma = 1, skew = skewFor, shape = shapeFor))
          var[i, ] <- -sigmaFor * (quantile)
          es[i, ] <- sigmaFor * ((shapeFor + quantile^2) / (shapeFor - 1) * (pdf_value / c))
        } else {
          stop("Unsupported distribution type. Use 'norm', 'std', or 'sstd'.")
        }
      } else {
        var[i, ] <- vol[i] * var_evt
        es[i, ] <- vol[i] * es_evt
      }
    } else {
      warning("xi, beta, or m / n_u invalid on date: ", as.character(df$Date[m + i - 1]))
    }
  }

  if (verbose > 0) {
    pb$terminate(); invisible()
  }

  # Create xts objects for VaR and ES
  dates <- tail(data$Date, n)
  VaR <- xts::xts(var, order.by = dates)
  ES <- xts::xts(es, order.by = dates)
  VOL <- xts::xts(vol, order.by = dates, colnames = "VOL")
  colnames(VOL) <- "VOL"
  results_xts <- xts::merge.xts(VaR, ES, VOL)

  if (any(is.na(results_xts))) {
    warning("There were NA values, carry forward last non-NA")
    results_xts <- zoo::na.locf(results_xts, na.rm = FALSE)
  }

  return(results_xts)
}

main <- function() {
  set.seed(123)
  df <- read.csv("C:/Users/chris/RStudioProjects/var-es-toolbox/data/clean_returns.csv", row.names = 1)
  asset_df <- df["DEBMc1"]
  asset_df$Date <- as.Date(rownames(asset_df))
  colnames(asset_df) <- c("Return", "Date")
  rownames(asset_df) <- NULL

  # Call the function
  result <- forecast_u_EVT_GARCH(
    asset_df,
    c = c(0.025, 0.05),
    n = 250,
    m = 250,
    p = 1,
    q = 1,
    r = 10,
    t = 0.9,
    verbose = 1,
    model = "sGARCH",
    dist = "std",
    solver = "hybrid",
    solver.control = list(trace = 0, tol = 1e-7)
  )

  print(tail(result))

}
