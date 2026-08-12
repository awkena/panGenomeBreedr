test_that("list_tables() works correctly in local mode", {
  # Explicitly test 'local' mode using the connection from setup.R
  tables <- list_tables(con = con_test, connect_db_mode = 'local')

  #  Check that the output is a character vector
  expect_type(tables, "character")

  # Check that it contains the expected tables 
  expect_true(all(c("annotations", "genotypes", "variants") %in% tables))
})

test_that("list_tables() defaults to 'local' mode correctly", {
  # Test that the default mode works as expected without explicit setting
  tables <- list_tables(con = con_test)
  expect_type(tables, "character")
  expect_true("variants" %in% tables)
})

test_that("list_tables() errors correctly in local mode without a valid connection", {
  # Test with a NULL connection
  expect_error(
    list_tables(con = NULL, connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode."
  )

  # Test with an explicitly closed connection
  closed_con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbDisconnect(closed_con, shutdown = TRUE)
  expect_error(
    list_tables(con = closed_con, connect_db_mode = 'local'),
    "A valid 'con' object is required for 'local' mode."
  )
})

test_that("list_tables() validates the connect_db_mode argument", {
  # Test with an invalid mode string, which should trigger an error from match.arg()
  expect_error(list_tables(con = con_test, connect_db_mode = "invalid_mode"))
})
