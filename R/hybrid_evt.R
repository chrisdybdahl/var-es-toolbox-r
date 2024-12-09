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

  spec <- rugarch::ugarchspec(
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
  for (i in 1:n) {
    if (verbose > 0) pb$tick()
    window_start <- i
    window_end <- m + i - 1
    window <- data_xts[window_start:window_end]

    fit <- rugarch::ugarchfit(
      spec = spec,
      data = window,
      solver = solver,
      solver.control = solver.control,
      ...
    )

    vol[i] <- rugarch::sigma(fit)[m]

    # Reparametrize for every r period
    if (i %% r == 0 || i == 1) {
      std_res <- rugarch::residuals(fit, standardize = TRUE)

      # Threshold selection
      u <- quantile(std_res, probs = t)

      # Calculate exceedances with a marked point process framework
      n_u <- length(std_res[std_res > u])

      if (n_u >= min_exceedances) {
        # Fit the Generalized Pareto Distribution with tryCatch
        evt_fit <- tryCatch({
          extRemes::fevd(std_res, threshold = u, type = "GP", method = "MLE")
        }, error = function(e) {
          warning("GPD fitting failed: ", e$message)
          return(NULL)
        })

        if (!any(is.null(evt_fit))) {
          # Extract parameters if fitting was successful
          xi <- evt_fit$results$par["shape"]
          beta <- evt_fit$results$par["scale"]
        } else {
          # Handle failed fitting
          warning("GPD fitting was not successful on date: ", as.character(df$Date[m + i - 1]), ". Retaining previous parameters.")
          xi <- xi
          beta <- beta
        }
      } else {
        warning("Not enough exceedances for EVT fitting on date: ", as.character(df$Date[m + i - 1]), ". Skipping EVT adjustment.")
        xi <- xi
        beta <- beta
      }
    }

    if (!any(is.na(c(xi, beta, m / n_u)))) {
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
      if (any(is.na(var_evt)) || any(is.na(es_evt))) {
        warning("VaR or ES calculation returned NA on date: ", as.character(df$Date[m + i - 1]))
      } else if (any(var_evt <= 0) || any(es_evt <= 0)) {
        warning("VaR or ES calculation returned non-positive values on date: ", as.character(df$Date[m + i - 1]))
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

