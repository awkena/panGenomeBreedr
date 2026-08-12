test_that("summarize_database() works correctly in local mode", {
  db_summary <- summarize_database(con = con_test, connect_db_mode = 'local')

  # Check output type and structure
  expect_s3_class(db_summary, "data.frame")
  expect_equal(colnames(db_summary), c("table", "n_rows"))

  # Check content
  expect_true(all(c("annotations", "genotypes", "metadata", "variants") %in% db_summary$table))
  expect_true(all(db_summary$n_rows >= 0))
  expect_true(is.numeric(db_summary$n_rows))
})

test_that("summarize_database() errors correctly without a valid connection", {
  expect_error(
    summarize_database(con = NULL, connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )
})

test_that("summarize_database() validates connect_db_mode argument", {
  expect_error(summarize_database(con = con_test, connect_db_mode = "invalid_mode"))
})
