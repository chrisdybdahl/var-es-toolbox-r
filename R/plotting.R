require(xts)
require(psych)
require(ggplot2)
require(gridExtra)
require(scales)
require(grid)

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

plot_var_es_predictions <- function(assets, conf_level, dir = "") {
  var_plots <- list()
  es_plots <- list()
  shared_legend <- NULL

  for (asset in assets) {
    # Construct file path and check for existence
    file_name <- file.path(paste0(dir, asset, "_predictions.csv"))
    if (!file.exists(file_name)) {
      warning("File not found: ", file_name)
      next
    }

    # Read CSV and parse dates
    predictions <- utils::read.csv(file_name, row.names = 1, stringsAsFactors = FALSE)

    # Attempt to parse row names as dates
    tryCatch({
      predictions$Date <- as.Date(rownames(predictions), format = "%d-%m-%y")
    }, error = function(e) {
      warning("Error parsing dates for asset ", asset, ": ", e$message)
      return(NULL)
    })

    print(predictions$Date)

    # Validate that Date column is now valid
    if (any(is.na(predictions$Date))) {
      warning("Invalid Date values detected for asset: ", asset)
      next
    }

    # Identify VaR and ES columns
    var_cols <- grep(paste0("\\.VaR\\.VaR_", conf_level, "$"), colnames(predictions), value = TRUE)
    es_cols <- grep(paste0("\\.ES\\.ES_", conf_level, "$"), colnames(predictions), value = TRUE)

    if (length(var_cols) == 0 || length(es_cols) == 0) {
      warning("No VaR/ES predictions found for asset: ", asset)
      next
    }

    # Extract model class names
    extract_model_class <- function(column_names) {
      gsub("_.*$", "", gsub("\\..*$", "", column_names))
    }
    var_models <- extract_model_class(var_cols)
    es_models <- extract_model_class(es_cols)

    # Create color palette for model classes
    unique_models <- unique(c(var_models, es_models))
    color_palette <- scales::hue_pal()(length(unique_models))
    model_colors <- setNames(color_palette, unique_models)

    # Prepare data for plotting
    var_data <- predictions[, c("Date", var_cols), drop = FALSE]
    es_data <- predictions[, c("Date", es_cols), drop = FALSE]

    var_long <- tidyr::pivot_longer(var_data, cols = -Date, names_to = "Model", values_to = "VaR")
    var_long$ModelName <- extract_model_class(var_long$Model)

    es_long <- tidyr::pivot_longer(es_data, cols = -Date, names_to = "Model", values_to = "ES")
    es_long$ModelName <- extract_model_class(es_long$Model)

    # Create individual VaR and ES plots
    var_plot <- ggplot2::ggplot(var_long, ggplot2::aes(x = Date, y = VaR, color = ModelName, group = Model)) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::scale_color_manual(values = model_colors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) # Improve date readability
      )

    es_plot <- ggplot2::ggplot(es_long, ggplot2::aes(x = Date, y = ES, color = ModelName, group = Model)) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::scale_color_manual(values = model_colors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    # Extract shared legend
    if (is.null(shared_legend)) {
      shared_legend <- ggplot2::ggplotGrob(var_plot + ggplot2::theme(legend.position = "bottom"))$grobs[[which(
        sapply(ggplot2::ggplotGrob(var_plot + ggplot2::theme(legend.position = "bottom"))$grobs, function(x) x$name == "guide-box"
        ))]]
    }

    # Store plots
    if (inherits(var_plot, "ggplot")) var_plots[[asset]] <- var_plot
    if (inherits(es_plot, "ggplot")) es_plots[[asset]] <- es_plot
  }

  # Filter valid plots
  var_plots <- Filter(function(x) inherits(x, "ggplot"), var_plots)
  es_plots <- Filter(function(x) inherits(x, "ggplot"), es_plots)

  # Combine plots and add legend
  if (length(var_plots) > 0) {
    grid.newpage()
    multiplot_var <- gridExtra::grid.arrange(
      gridExtra::arrangeGrob(grobs = var_plots, ncol = 2),
      shared_legend,
      ncol = 1,
      heights = c(9, 1)
    )
    print(multiplot_var)
  }

  if (length(es_plots) > 0) {
    grid.newpage()
    multiplot_es <- gridExtra::grid.arrange(
      gridExtra::arrangeGrob(grobs = es_plots, ncol = 2),
      shared_legend,
      ncol = 1,
      heights = c(9, 1)
    )
    print(multiplot_es)
  }
}
