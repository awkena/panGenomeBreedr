test_that("list_sqlite_tables() correctly lists tables in the database", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Use the function to list tables from your active mini Parquet views
  tables <- list_sqlite_tables(con_test)

  # Check structure
  expect_type(tables, "character")

  # Expect all 4 essential tables to be mounted and discovered
  expect_true(all(
    c("genotypes", "variants", "metadata", "annotations") %in% tables
  ))
  expect_length(tables, 4)
})