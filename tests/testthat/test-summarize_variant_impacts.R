test_that("summarize_variant_impacts() works in local mode", {
  impact_summary <- summarize_variant_impacts(con = con_test, connect_db_mode = 'local')

  # Check that the output is a data frame and has rows
  expect_s3_class(impact_summary, "data.frame")
  expect_gt(nrow(impact_summary), 0)

  # Check for expected column structure
  expect_true("chrom" %in% colnames(impact_summary))
  expect_true(any(grepl("^impact_", colnames(impact_summary))))

  # Check that NAs were correctly replaced with 0
  numeric_cols <- impact_summary[, -which(colnames(impact_summary) == "chrom")]
  expect_false(any(is.na(numeric_cols)))
})

test_that("summarize_variant_impacts() errors correctly without a valid connection", {
  expect_error(
    summarize_variant_impacts(con = NULL, connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("summarize_variant_impacts() validates connect_db_mode argument", {
  expect_error(summarize_variant_impacts(con = con_test, connect_db_mode = "invalid_mode"))
})
