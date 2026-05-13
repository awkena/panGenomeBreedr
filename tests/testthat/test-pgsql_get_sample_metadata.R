test_that("pgsql_get_sample_metadata logic holds in SQLite (Mock-independent)", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # Ensure connection is closed cleanly
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Create a dummy metadata table that includes the array_index column
  DBI::dbExecute(
    con,
    "CREATE TABLE metadata (lib TEXT, countryorigin TEXT, array_index INTEGER)"
  )

  # Insert mock data with the array index
  DBI::dbExecute(
    con,
    "INSERT INTO metadata VALUES ('Sample1', 'Ghana', 1), ('Sample2', 'USA', 2)"
  )

  # Test filtering logic
  res <- pgsql_get_sample_metadata(
    con,
    query_col = "countryorigin",
    query_value = "Ghana"
  )

  expect_equal(nrow(res), 1)
  expect_equal(res$lib, "Sample1")
})