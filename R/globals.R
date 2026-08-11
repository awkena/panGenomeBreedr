
# Declare global variables to prevent R CMD check NOTEs caused by ggplot2 NSE (Non-Standard Evaluation)
utils::globalVariables(c(
  "Start", "End", "ymin", "ymax", "Feature",
  "x_start", "x_end", "plot_y", "impact",
  "variant_type", "pos.x", "clean_label"
))
