require(esback)
require(Rcpp)
require(rugarch)
require(segMGarch)
require(Rcpp)

# Kupiec (1995) and Christoffersen (1998)
run_kc_backtest <- function(actual, VaR, alpha) {
  result <- rugarch::VaRTest(alpha = alpha, actual = actual, VaR = VaR)
  return(list(
    ExpectedExceedances = result$expected.exceed,
    ActualExceedances = result$actual.exceed,
    UC = result$uc.LRp,
    CC = result$cc.LRp
  ))
}

# Christoffersen and Pelletier (2004)
# TODO: Remember the following note
# Following Christoffersen and Pelletier (2004), the Weibull distribution is used with parameter ‘b=1’ representing the
# case of the exponential. A future release will include the choice of using a bootstrap method to evaluate the p-value,
# and until then care should be taken when evaluating series of length less than 1000 as a rule of thumb.
run_vdt_backtest <- function(actual, VaR, alpha) {
  result <- rugarch::VaRDurTest(alpha = alpha, actual = actual, VaR = VaR)
  return(result$LRp)
}

# Dynamic Quantile Test (Engle and Manganelli, 2004)
run_dq_backtest <- function(actual, VaR, alpha, lag = 1, lag_hit = 1, lag_var = 1) {
  result <- segMGarch::DQtest(y = as.numeric(actual), VaR = VaR, VaR_level = alpha, lag = lag, lag_hit = lag_hit, lag_var = lag_var)
  return(result[[2]])
}

# McNeil and Frey (2000) test
run_er_backtest <- function(actual, VaR, ES, vol = NULL, b = 1000) {
  result <- esback::er_backtest(r = actual, q = VaR, e = ES, s = vol, B = b)
  return(if (is.null(vol)) result$pvalue_twosided_simple else result$pvalue_twosided_standardized)
}

# McNeil and Frey (2000) test second version
run_er_backtest_2 <- function(actual, VaR, ES, alpha, boot = FALSE, b = 1000) {
  result <- rugarch::ESTest(alpha = alpha, actual = actual, VaR = VaR, ES = ES, boot = boot, n.boot = b)
  return(result$p.value)
}

# Nolde and Ziegel (2017) test
run_coc_backtest <- function(actual, VaR, ES, alpha, vol = NULL) {
  result <- cc_backtest(r = actual, q = VaR, e = ES, s = vol, alpha = alpha)
    return(if (is.null(vol)) result$pvalue_twosided_simple else result$pvalue_twosided_general)
}

# Bayer and Dimitriadis (2020) tests
run_esr_backtests <- function(actual, VaR, ES, alpha, b = 0) {
  # Perform ESR backtests for all three versions
  BD_results <- lapply(1:3, function(version) {
    tryCatch(
      {
        BD <- esback::esr_backtest(r = actual, q = VaR, e = ES, alpha = alpha, version = version, B = b)
        return(BD$pvalue_twosided_asymptotic)
      },
      error = function(e) NA
    )
  })

  # Return as a named list for each version
  return(list(
    ESR1 = BD_results[[1]],
    ESR2 = BD_results[[2]],
    ESR3 = BD_results[[3]]
  ))
}

run_backtests <- function(actual, VaR, ES, alpha, prefix, VOL = NULL, b = 1000) {
  # Initialize a data frame with NA values for all expected outputs
  backtest_template <- data.frame(
    ExpectedExceedances = NA,
    ActualExceedances = NA,
    UC = NA,
    CC = NA,
    VDT = NA,
    DQ = NA,
    ER_Simple = NA,
    ER_Standardized = NA,
    ER_2 = NA,
    CoC_Simple = NA,
    CoC_Generalized = NA,
    ESR1 = NA,
    ESR2 = NA,
    ESR3 = NA
  )

  # Use tryCatch to handle errors gracefully
  safe_run <- function(expr) {
    tryCatch(
      expr,
      error = function(e) {
        message("Error encountered: ", e$message)
        return(NA)
      }
    )
  }

  # Kupiec (1995) and Christoffersen (1998) backtest
  kc_results <- safe_run(run_kc_backtest(actual = actual, VaR = VaR, alpha = alpha))
  backtest_template$ExpectedExceedances <- ifelse(is.list(kc_results), kc_results$ExpectedExceedances, NA)
  backtest_template$ActualExceedances <- ifelse(is.list(kc_results), kc_results$ActualExceedances, NA)
  backtest_template$UC <- ifelse(is.list(kc_results), kc_results$UC, NA)
  backtest_template$CC <- ifelse(is.list(kc_results), kc_results$CC, NA)

  # Christoffersen and Pelletier (2004) backtest
  vdt_pvalue <- safe_run(run_vdt_backtest(actual = actual, VaR = VaR, alpha = alpha))
  backtest_template$VDT <- vdt_pvalue

  # Dynamic Quantile (DQ) test
  dq_pvalue <- safe_run(run_dq_backtest(actual = actual, VaR = VaR, alpha = alpha))
  backtest_template$DQ <- dq_pvalue

  # McNeil and Frey (2000) backtest (simple and standardized)
  er_simple <- safe_run(run_er_backtest(actual = actual, VaR = VaR, ES = ES, vol = NULL, b = b))
  er_standardized <- if (!is.null(VOL)) {
    safe_run(run_er_backtest(actual = actual, VaR = VaR, ES = ES, vol = VOL, b = b))
  } else {
    "Not applicable"
  }
  backtest_template$ER_Simple <- er_simple
  backtest_template$ER_Standardized <- er_standardized

  # McNeil and Frey (2000) test second version
  er_pvalue_2 <- safe_run(run_er_backtest_2(actual = actual, VaR = VaR, ES = ES, alpha = alpha))
  backtest_template$ER_2 <- er_pvalue_2

  # Nolde and Ziegel (2017) backtest (simple and standardized)
  coc_simple <- safe_run(run_coc_backtest(actual = actual, VaR = VaR, ES = ES, alpha = alpha, vol = NULL))
  coc_generalized <- if (!is.null(VOL)) {
    safe_run(run_coc_backtest(actual = actual, VaR = VaR, ES = ES, alpha = alpha, vol = VOL))
  } else {
    "Not applicable"
  }
  backtest_template$CoC_Simple <- coc_simple
  backtest_template$CoC_Generalized <- coc_generalized

  # Bayer and Dimitriadis (2020) backtests
  esr_results <- safe_run(run_esr_backtests(actual = actual, VaR = VaR, ES = ES, alpha = alpha))
  backtest_template$ESR1 <- ifelse(is.list(esr_results), esr_results$ESR1, NA)
  backtest_template$ESR2 <- ifelse(is.list(esr_results), esr_results$ESR2, NA)
  backtest_template$ESR3 <- ifelse(is.list(esr_results), esr_results$ESR3, NA)

  row.names(backtest_template) <- prefix

  return(backtest_template)
}

# Negative Asymmetric Laplace Log-Likelihood (AL) Loss
custom_loss <- function(returns, VaR, ES, c) {
  loss_vector <- -log((c - 1) / ES) - ((returns - VaR) * (c - (returns <= VaR))) / (c * ES)
  return(loss_vector)
}

cppFunction('
NumericVector custom_loss_rcpp(NumericVector returns, NumericVector VaR, NumericVector ES, double c) {
  int n = returns.size();
  NumericVector loss_vector(n);

  for (int i = 0; i < n; i++) {
    if (ES[i] == 0) {
      loss_vector[i] = NA_REAL; // Handle cases where ES is zero to avoid division by zero
    } else {
      double term1 = -log((c - 1) / ES[i]);
      double term2 = ((returns[i] - VaR[i]) * (c - (returns[i] <= VaR[i]))) / (c * ES[i]);
      loss_vector[i] = term1 - term2;
    }
  }

  return loss_vector;
}
')

compute_losses_and_backtests <- function(df, models, c_levels, n, m, b = 1000, asset = "results") {
  # File paths for saving results
  predictions_file <- paste0(asset, "_predictions.csv")
  losses_file <- paste0(asset, "_losses.csv")
  backtests_file <- paste0(asset, "_backtests.csv")

  # Load existing predictions if available
  if (file.exists(predictions_file)) {
    existing_predictions <- read.csv(predictions_file)
    print(existing_predictions)
  } else {
    existing_predictions <- data.frame(Date = tail(df$Date, n))
  }

  # Load existing losses and backtests if available
  if (FALSE && file.exists(losses_file)) {
    losses_df <- read.csv(losses_file, row.names = 1)
  } else {
    losses_df <- data.frame(Date = tail(df$Date, n))
  }

  if (FALSE && file.exists(backtests_file)) {
    backtests_df <- read.csv(backtests_file, row.names = 1)
  } else {
    backtests_df <- data.frame()
  }

  # Identify models with predictions
  models_with_predictions <- if (ncol(existing_predictions) > 1) {
    unique(sub("\\..*", "", colnames(existing_predictions)))
  } else {
    character(0)
  }

  # Initialize backtest template
  backtest_template <- data.frame(
    ExpectedExceedances = NA,
    ActualExceedances = NA,
    UC = NA,
    CC = NA,
    VDT = NA,
    DQ = NA,
    ER_Simple = NA,
    ER_Standardized = NA,
    ER_2 = NA,
    CoC_Simple = NA,
    CoC_Generalized = NA,
    ESR1 = NA,
    ESR2 = NA,
    ESR3 = NA
  )

  # Loop through each model
  for (model_name in names(models)) {
    print(paste("Processing model:", model_name))
    compute <- TRUE
    successfull <- TRUE

    # Check if predictions for this model exist
    if (model_name %in% models_with_predictions) {
      # Reuse predictions from the existing file
      model_predictions <- existing_predictions[, grep(paste0("^", model_name, "\\."), colnames(existing_predictions)), drop = FALSE]
      message(paste("Using existing predictions for model:", model_name))

      # Ensure both VaR and ES predictions are available
      var_available <- any(grepl("VaR", colnames(model_predictions)))
      es_available <- any(grepl("ES", colnames(model_predictions)))
      if (!var_available || !es_available) {
        message(paste("Incomplete predictions for model:", model_name, "- Recomputing"))
      } else {
        compute <- FALSE
      }
    }

    if (compute) {
      # Compute new predictions if not already available
      model_forecasts <- tryCatch(models[[model_name]](df, c_levels, n, m),
        error = function(e) {
          message(paste("Model", model_name, "failed:", e$message))
          return(NULL)
        }
      )

      if (!is.null(model_forecasts)) {
        # Ensure the forecasts are in matrix format
        model_forecasts <- as.matrix(model_forecasts)
        var_cols <- grep("^VaR_", colnames(model_forecasts))
        es_cols <- grep("^ES_", colnames(model_forecasts))

        # Ensure both VaR and ES columns are available
        if (length(var_cols) != length(c_levels) || length(es_cols) != length(c_levels)) {
          message(paste("Incomplete VaR/ES columns for model:", model_name, "- Skipping"))
          next
        }

        for (c_idx in seq_along(c_levels)) {
          c <- c_levels[c_idx]
          existing_predictions[, paste0(model_name, ".VaR.VaR_", c)] <- model_forecasts[, var_cols[c_idx]]
          existing_predictions[, paste0(model_name, ".ES.ES_", c)] <- model_forecasts[, es_cols[c_idx]]
        }

      } else {
        successfull <- FALSE
      }
    }

    # Compute losses and backtests using the predictions
    for (c_idx in seq_along(c_levels)) {
      c <- c_levels[c_idx]

      if (successfull) {
        VaR_col <- existing_predictions[, paste0(model_name, ".VaR.VaR_", c), drop = TRUE]
        ES_col <- existing_predictions[, paste0(model_name, ".ES.ES_", c), drop = TRUE]

        # Compute losses
        loss <- tryCatch(
          custom_loss_rcpp(tail(df$Return, n), -VaR_col,-ES_col, c),
          error = function(e) {
            message(paste("Loss computation failed for model:", model_name, "confidence level:", c))
            return(rep(NA, n))
          }
        )

        # Run backtests and add average loss
        backtest <- tryCatch(
          {
            backtest <- run_backtests(
              actual = tail(df$Return, n),
              VaR = -VaR_col,
              ES = -ES_col,
              alpha = c,
              prefix = paste0(model_name, "_", c),
              b = b
            )
            backtest
          },
          error = function(e) {
            message(paste("Backtest failed for model:", model_name, "confidence level:", c))
            backtest_template
          }
        )
      } else {
        loss <- rep(NA, n)
        backtest <- backtest_template
      }

      # Store losses
      losses_df[[paste0(model_name, "_", c)]] <- loss

      # Compute average loss
      avg_loss <- mean(loss, na.rm = TRUE)

      backtest$AverageLoss <- avg_loss
      rownames(backtest) <- paste0(model_name, "_", c)

      # Store backtests
      backtests_df <- rbind(backtests_df, backtest)
    }

    # Save losses, backtests, and predictions continuously
    write.csv(losses_df, losses_file, row.names = FALSE)
    write.csv(backtests_df, backtests_file, row.names = TRUE)
    write.csv(existing_predictions, predictions_file, row.names = FALSE)
  }

  return(list(
    predictions = existing_predictions,
    losses = losses_df,
    backtests = backtests_df
  ))
}


# Function to perform MCS test for each confidence level
perform_mcs_test_by_confidence <- function(losses, alpha_levels, nboot) {
  mcs_results <- list()

  # Exclude the "Date" column and process only model-related columns
  relevant_columns <- colnames(losses)[colnames(losses) != "Date"]

  # Ensure there are relevant columns for processing
  if (length(relevant_columns) == 0) {
    warning("No relevant columns found for MCS test.")
    return(mcs_results)
  }

  for (conf_level in unique(sub(".*_", "", relevant_columns))) {
    print(paste("Processing confidence level:", conf_level))

    # Filter losses for current confidence level
    conf_losses <- losses[, grepl(paste0("_", conf_level, "$"), relevant_columns), drop = FALSE]

    if (ncol(conf_losses) == 0) {
      warning(paste("No columns found for confidence level:", conf_level))
      next
    }

    # Replace NA and run MCS test
    conf_losses[is.na(conf_losses)] <- 1e6

    for (alpha in alpha_levels) {
      tryCatch({
        mcs_result <- rugarch::mcsTest(conf_losses, alpha = alpha, nboot = nboot, boot = "stationary")

        includedR_indices <- mcs_result$includedR
        includedSQ_indices <- mcs_result$includedSQ

        mcs_results[[paste0("alpha_", alpha, "_conf_", conf_level)]] <- list(
          includedR_names = colnames(conf_losses)[includedR_indices],
          excludedR_names = colnames(conf_losses)[setdiff(seq_len(ncol(conf_losses)), includedR_indices)],
          includedSQ_names = colnames(conf_losses)[includedSQ_indices],
          excludedSQ_names = colnames(conf_losses)[setdiff(seq_len(ncol(conf_losses)), includedSQ_indices)]
        )
      }, error = function(e) {
        message(paste("MCS test failed for alpha:", alpha, "confidence level:", conf_level, "-", e$message))
      })
    }
  }

  return(mcs_results)

}
