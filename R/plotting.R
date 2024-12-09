require(xts)
require(rugarch)
require(psych)
require(ggplot2)
require(gridExtra)
require(tseries)

# Function to plot all returns
plot_all_returns <- function(data) {
  # Assume the relevant columns are identified (excluding "Date")
  return_cols <- setdiff(names(data), "Date")

  # Calculate global y-axis limits
  y_limits <- range(data[return_cols], na.rm = TRUE)

  # Create a list to store individual plots
  plot_list <- lapply(return_cols, function(col) {
    ggplot(data, aes_string(x = "Date", y = col)) +
      geom_line(color = "black") +
      labs(title = col, x = "", y = "") +
      theme_minimal() +
      scale_y_continuous(limits = y_limits) # Apply consistent y-axis limits
  })

  # Combine all plots into a grid
  do.call(grid.arrange, c(plot_list, ncol = 2)) # Adjust ncol to control layout
}

# Function to create QQ plots for all columns in a dataframe
plot_all_qq <- function(data) {
  # Exclude 'Date' column
  data_cols <- setdiff(names(data), "Date")

  # Calculate global y-axis limits for QQ plots
  qq_limits <- quantile(unlist(data[data_cols]), probs = c(0.01, 0.99), na.rm = TRUE)

  # Create a list to store individual QQ plots
  plot_list <- lapply(data_cols, function(col) {
    ggplot(data, aes_string(sample = col)) +
      stat_qq() +
      stat_qq_line() +
      labs(title = col, x = "", y = "") +
      theme_minimal() +
      scale_y_continuous(limits = qq_limits) # Apply consistent y-axis limits
  })

  # Combine all QQ plots into a grid
  do.call(grid.arrange, c(plot_list, ncol = 2)) # Adjust ncol to control layout
}

# Function to plot all prices
plot_all_timeseries <- function(data) {
  timeseries_cols <- setdiff(names(data), "Date")

  # Define a custom matte color palette
  custom_colors <- c("#597BAF", "#E6A157", "#7F9C66", "#C7615D",
                     "#B39EB5", "#A07C72", "#D69AB3", "#A6A6A6")

  # Map colors to column names
  color_mapping <- custom_colors[seq_along(timeseries_cols)]
  names(color_mapping) <- timeseries_cols

  # Create the base plot
  p <- ggplot(data, aes(x = Date)) +
    labs(x = "Date", y = "Value") +
    theme_minimal()

  # Add each time series as a separate line
  for (col in timeseries_cols) {
    p <- p + geom_line(aes_string(y = col, color = shQuote(col)))
  }

  # Finalize the plot with consistent y-axis limits
  p + scale_color_manual(values = color_mapping, name = "Series")
}

plot_all_timeseries <- function(data) {
  timeseries_cols <- setdiff(names(data), "Date")

  # Define a custom matte color palette
  custom_colors <- c("#597BAF", "#E6A157", "#7F9C66", "#C7615D",
                     "#B39EB5", "#A07C72", "#D69AB3", "#A6A6A6")

  # Map colors to column names
  color_mapping <- custom_colors[seq_along(timeseries_cols)]
  names(color_mapping) <- timeseries_cols

  # Create the base plot
  p <- ggplot(data, aes(x = Date)) +
    labs(x = "Date", y = "Value") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove default legend

  # Add each time series as a separate line
  for (col in timeseries_cols) {
    p <- p + geom_line(aes_string(y = col, color = shQuote(col)))
  }

  # Prepare legend positions and labels
  legend_df <- data.frame(
    x = rep(max(data$Date) - (max(data$Date) - min(data$Date)) * 0.1, length(timeseries_cols)),
    y = seq(
      from = max(unlist(data[timeseries_cols], use.names = FALSE)) * 0.9,
      to = max(unlist(data[timeseries_cols], use.names = FALSE)) * 0.7,
      length.out = length(timeseries_cols)
    ),
    label = names(color_mapping),
    color = color_mapping
  )

  # Add custom legend inside the plot
  p <- p +
    geom_text(data = legend_df, aes(x = x, y = y, label = label, color = label), hjust = 0, size = 3) +
    geom_point(data = legend_df, aes(x = x - (max(data$Date) - min(data$Date)) * 0.01, y = y, color = label), size = 3)

  # Finalize the plot with custom colors
  p + scale_color_manual(values = color_mapping, name = "Series")
}

# Function to plot base vs. peak prices in subplots
plot_base_vs_peak <- function(data, base_cols, peak_cols) {
  if (length(base_cols) != length(peak_cols)) {
    stop("Base and peak column lists must have the same length.")
  }

  # Calculate global y-axis limits
  all_cols <- c(base_cols, peak_cols)
  y_limits <- range(data[all_cols], na.rm = TRUE)

  # Create a list to store individual plots
  plot_list <- lapply(seq_along(base_cols), function(i) {
    base_col <- base_cols[i]
    peak_col <- peak_cols[i]

    ggplot(data, aes(x = Date)) +
      geom_line(aes_string(y = base_col), color = "#597BAF", linetype = "solid") +
      geom_line(aes_string(y = peak_col), color = "#E6A157", linetype = "dashed") +
      labs(title = paste(base_col, "vs", peak_col),
           x = "",
           y = "") +
      theme_minimal() +
      scale_y_continuous(limits = y_limits) # Apply consistent y-axis limits
  })

  # Combine all plots into a grid
  do.call(grid.arrange, c(plot_list, ncol = 2)) # Adjust ncol to control layout
}

# Define the function for plotting returns, VaR, and ES
plot_var_es <- function(dates, returns, var, es) {
  # Create a data frame for plotting
  plot_data <- data.frame(
    Date = dates,
    Return = returns,
    VaR = var,
    ES = es
  )

  # Generate the plot using ggplot2 with labels
  ggplot(plot_data, aes(x = Date)) +
    geom_line(aes(y = Return, color = "Return"), linewidth = 0.7, alpha = 0.8) +    # Actual returns
    geom_line(aes(y = VaR, color = "VaR"), linetype = "dashed", linewidth = 0.7) + # VaR
    geom_line(aes(y = ES, color = "ES"), linetype = "dotted", linewidth = 0.7) +   # ES
    ggtitle("Actual Returns vs. VaR and ES") +
    xlab("Date") +
    ylab("Returns / VaR / ES") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    scale_color_manual(values = c("Return" = "blue", "VaR" = "red", "ES" = "darkorange")) +
    scale_y_continuous(labels = scales::percent) +
    labs(color = "Legend")
}

plot_var_es_predictions <- function(assets, conf_level) {
  # Initialize lists to store plots
  var_plots <- list()
  es_plots <- list()

  for (asset in assets) {
    # Read predictions CSV
    file_name <- file.path("results", paste0(asset, "_predictions.csv"))
    if (!file.exists(file_name)) {
      warning(paste("File not found:", file_name))
      next
    }
    predictions <- utils::read.csv(file_name, row.names = 1, stringsAsFactors = FALSE)
    predictions$Date <- as.Date(rownames(predictions))

    # Extract relevant VaR and ES columns
    var_cols <- grep(paste0("\\.VaR\\.VaR_", conf_level, "$"), colnames(predictions), value = TRUE)
    es_cols <- grep(paste0("\\.ES\\.ES_", conf_level, "$"), colnames(predictions), value = TRUE)

    if (length(var_cols) == 0 || length(es_cols) == 0) {
      warning(paste("No VaR/ES predictions found for confidence level", conf_level, "in", asset))
      next
    }

    # Reshape data for ggplot using tidyr
    var_data <- predictions[, c("Date", var_cols), drop = FALSE]
    es_data <- predictions[, c("Date", es_cols), drop = FALSE]

    var_long <- tidyr::pivot_longer(var_data, cols = -Date, names_to = "Model", values_to = "VaR")
    es_long <- tidyr::pivot_longer(es_data, cols = -Date, names_to = "Model", values_to = "ES")

    # Plot VaR without labels
    var_plot <- ggplot2::ggplot(var_long, ggplot2::aes(x = Date, y = VaR, group = Model)) +
      ggplot2::geom_line(alpha = 0.5) +
      ggplot2::labs(title = paste(asset, "VaR Predictions (", conf_level, ")"), x = "Date", y = "VaR") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")  # Remove legend

    # Plot ES without labels
    es_plot <- ggplot2::ggplot(es_long, ggplot2::aes(x = Date, y = ES, group = Model)) +
      ggplot2::geom_line(alpha = 0.5) +
      ggplot2::labs(title = paste(asset, "ES Predictions (", conf_level, ")"), x = "Date", y = "ES") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")  # Remove legend

    # Append plots to the list
    var_plots[[asset]] <- var_plot
    es_plots[[asset]] <- es_plot
  }

  # Arrange and display plots
  multiplot_var <- gridExtra::grid.arrange(grobs = var_plots, ncol = 1, top = "VaR Predictions")
  multiplot_es <- gridExtra::grid.arrange(grobs = es_plots, ncol = 1, top = "ES Predictions")

  print(multiplot_var)
  print(multiplot_es)
}
