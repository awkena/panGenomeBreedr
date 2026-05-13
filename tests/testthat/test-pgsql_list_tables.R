test_that("pgsql_list_tables returns expected table names from mock DB", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")

  # Explicitly point to the fixtures folder
  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      # Connect to database
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # List tables
      table_names <- pgsql_list_tables(con)

      # Validate results
      expect_type(table_names, "character")
      expect_true("variants" %in% table_names)
      expect_true("metadata" %in% table_names)

      # Disconnect
      DBI::dbDisconnect(con)
    })
  })
})

test_that("pgsql_list_tables throws error on invalid connection", {
  # Create connection and close to invalidate
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbDisconnect(con)

  # Check error handling
  expect_error(
    pgsql_list_tables(con),
    "The provided database connection is not valid"
  )
})