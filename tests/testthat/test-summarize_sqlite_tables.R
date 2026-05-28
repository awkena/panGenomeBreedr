test_that("summarize_sqlite_tables() returns correct table names and row counts", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Run the function against your active mini Parquet tables
  summary_df <- summarize_sqlite_tables(con_test)

  # Validate output structure
  expect_true(is.data.frame(summary_df))
  expect_named(summary_df, c("table", "n_rows"))
  expect_equal(nrow(summary_df), 4) # Should discover genotypes, variants, metadata, and annotations

  # Check that the discovered table labels are valid schema members
  expect_true(all(
    c("variants", "annotations", "genotypes", "metadata") %in% summary_df$table
  ))

  # Assert that reported row counts are non-negative numeric vectors
  expect_type(summary_df$n_rows, "double") # DuckDB count queries return numeric types
  expect_true(all(summary_df$n_rows >= 0))
})