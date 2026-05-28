test_that("variant_impact_summary() returns impact summary in wide format", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Run the function against your active mini Parquet dataset
  summary_df <- variant_impact_summary(con_test)

  # Check basic structure properties
  expect_true(is.data.frame(summary_df))
  expect_true("chrom" %in% names(summary_df))

  # Verify that wide-format reshaping prefixing works flawlessly
  impact_cols <- names(summary_df)[grepl("^impact_", names(summary_df))]
  expect_true(length(impact_cols) > 0)

  # Check that it returns structured data lines for chromosomes present in your slice
  expect_true(nrow(summary_df) >= 1)
  expect_true("Chr05" %in% summary_df$chrom)

  # Ensure the cell values inside pivoted columns are numeric count vectors
  # and do not contain raw unhandled NA voids
  for (col in impact_cols) {
    expect_type(summary_df[[col]], "double")
    expect_true(all(!is.na(summary_df[[col]])))
    expect_true(all(summary_df[[col]] >= 0))
  }
})