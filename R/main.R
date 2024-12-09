# Set directories
data_dir <- file.path("data")
backtest_path <- file.path("R", "backtesting.R")
analysis_path <- file.path("R", "analysis.R")
plotting_path <- file.path("R", "plotting.R")
hs_path <- file.path("R", "hs.R")
garch_path <- file.path("R", "garch.R")
hybrid_evt_path <- file.path("R", "hybrid_evt.R")
gas_path <- file.path("R", "gas.R")
caviar_path <- file.path("R", "caviar.R")

source(backtest_path)
source(analysis_path)
source(plotting_path)
source(hs_path)
source(garch_path)
source(hybrid_evt_path)
source(gas_path)
source(caviar_path)

df <- read.csv("data/clean_returns.csv", row.names = 1)

# Define example models
models <- list(
  "HS" = function(df, c, n, m) forecast_u_HS(df, c = c, n = n, m = m),
  "FHS_EWMA" = function(df, c, n, m) forecast_u_FHS_EWMA(df, c = c, n = n, m = m),
  "FHS_sGARCH_norm" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "norm"),
  "FHS_sGARCH_std" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "std"),
  "FHS_sGARCH_sstd" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "sstd"),
  "FHS_gjrGARCH_norm" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "norm"),
  "FHS_gjrGARCH_std" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "std"),
  "FHS_gjrGARCH_sstd" = function(df, c, n, m) RollFHS(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "sstd"),
  "sGARCH_norm" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "norm"),
  "sGARCH_std" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "std"),
  "sGARCH_sstd" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "sGARCH", dist = "sstd"),
  "gjrGARCH_norm" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "norm"),
  "gjrGARCH_std" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "std"),
  "gjrGARCH_sstd" = function(df, c, n, m) forecast_u_GARCH(df, c = c, n = n, m = m, r = 10, model = "gjrGARCH", dist = "sstd"),
  "EVT_sGARCH_norm" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "sGARCH", dist = "norm"),
  "EVT_sGARCH_std" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "sGARCH", dist = "std"),
  "EVT_sGARCH_sstd" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "sGARCH", dist = "sstd"),
  "EVT_gjrGARCH_norm" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "gjrGARCH", dist = "norm"),
  "EVT_gjrGARCH_std" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "gjrGARCH", dist = "std"),
  "EVT_gjrGARCH_sstd" = function(df, c, n, m) forecast_u_EVT_GARCH(df, c = c, n = n, m = m, t = 0.9, r = 10, model = "gjrGARCH", dist = "sstd"),
  "GAS_norm" = function(df,c, n, m) forecast_u_GAS(df, c = c, n = n, m = m, r = 10, dist = "norm"),
  "GAS_std" = function(df,c, n, m) forecast_u_GAS(df, c = c, n = n, m = m, r = 10, dist = "std"),
  "GAS_sstd" = function(df,c, n, m) forecast_u_GAS(df, c = c, n = n, m = m, r = 10, dist = "sstd"),
  "CAViaR_ADAPTIVE_AR" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "ADAPTIVE", es_model = "AR"),
  "CAViaR_SAV_AR" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "SAV", es_model = "AR"),
  "CAViaR_AS_AR" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "AS", es_model = "AR"),
  "CAViaR_indirectGARCH_AR" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "indirectGARCH", es_model = "AR"),
  "CAViaR_ADAPTIVE_MULT" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "ADAPTIVE", es_model = "MULT"),
  "CAViaR_SAV_MULT" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "SAV", es_model = "MULT"),
  "CAViaR_AS_MULT" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "AS", es_model = "MULT"),
  "CAViaR_indirectGARCH_MULT" = function(df, c, n, m) forecast_u_CAViaR(df, c = c, n = n, m = m, r = 10, var_model = "indirectGARCH", es_model = "MULT")
)


# Set parameters
c <- c(0.025, 0.05)
n <- 989
m <- 250
p <- 1
q <- 1
t <- 0.9
refit <- 10
alpha_levels <- c(0.1, 0.05)
nboot <- 1000

assets <- setdiff(colnames(return_df), "Date")
#assets <- c("DEBMc1")

options(warn = 1)

for (asset in assets) {
  asset_df <- df[asset]
  asset_df$Date <- as.Date(rownames(asset_df))
  colnames(asset_df) <- c("Return", "Date")
  rownames(asset_df) <- NULL

  print(tail(asset_df,10))
  losses_and_backtests <- compute_losses_and_backtests(asset_df, models, c, n, m, asset = asset)

  losses <- losses_and_backtests$losses

  mcs_results_by_confidence <- perform_mcs_test_by_confidence(losses, alpha_levels, nboot)

  # Organize results into a data frame for easier interpretation
  mcs_table <- do.call(rbind, lapply(names(mcs_results_by_confidence), function(name) {
    result <- mcs_results_by_confidence[[name]]
    data.frame(
      Alpha_Confidence = name,
      IncludedR = paste(result$includedR_names, collapse = ", "),
      ExcludedR = paste(result$excludedR_names, collapse = ", "),
      IncludedSQ = paste(result$includedSQ_names, collapse = ", "),
      ExcludedSQ = paste(result$excludedSQ_names, collapse = ", ")
    )
  }))

  write.csv(mcs_table, paste0(asset, "_mcs.csv"), row.names = FALSE)

}

plots <- plot_var_es_predictions("DBc1", 0.05, dir = "")

