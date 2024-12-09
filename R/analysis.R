require(psych)
require(dplyr)


compute_yearly_std <- function(df) {
  # Ensure the Date column is in Date format
  if (!"Date" %in% colnames(df)) {
    stop("The CSV file must have a column named 'Date'.")
  }
  df$Date <- as.Date(df$Date)

  # Extract the year from the Date column
  df$Year <- format(df$Date, "%Y")

  # Group by year and compute standard deviation for each column
  # Exclude the Date and Year columns from calculations
  numeric_columns <- setdiff(names(df), c("Date", "Year"))
  yearly_std <- dplyr::group_by(df, Year) %>%
    dplyr::summarise(across(all_of(numeric_columns), ~ psych::sd(.x, na.rm = TRUE)))

  # Convert back to data frame and return
  yearly_std <- as.data.frame(yearly_std)

  return(yearly_std)
}
