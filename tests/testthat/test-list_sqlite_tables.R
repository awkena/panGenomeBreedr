test_that("list_sqlite_tables() correctly lists tables in a SQLite database", {
  # Skip test on CRAN and in case required packages aren't available
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create temporary SQLite database
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Add two example tables
  DBI::dbWriteTable(con, "genes", data.frame(id = 1:3, name = c("A", "B", "C")))
  DBI::dbWriteTable(con, "variants", data.frame(pos = c(100, 200), allele = c("A", "T")))

  DBI::dbDisconnect(con)

  # Use your function to list tables
  tables <- list_sqlite_tables(db_path)

  # Expect both tables to be found
  expect_true(all(c("genes", "variants") %in% tables))
  expect_length(tables, 2)

  # Cleanup
  unlink(db_path)
})
