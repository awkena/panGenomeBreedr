test_that("summarize_sqlite_tables() returns correct table names and row counts", {
  skip_on_cran()
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  # Create a temporary SQLite database
  db_path <- tempfile(fileext = ".db")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Create multiple tables with known row counts
  DBI::dbWriteTable(con, "variants", data.frame(id = 1:3, value = letters[1:3]))  # 3 rows
  DBI::dbWriteTable(con, "annotations", data.frame(id = 1:2, type = c("MODERATE", "HIGH")))  # 2 rows
  DBI::dbWriteTable(con, "meta", data.frame(info = character(0)))  # 0 rows

  DBI::dbDisconnect(con)

  # Run the function
  summary_df <- summarize_sqlite_tables(db_path)

  # Validate output
  expect_true(is.data.frame(summary_df))
  expect_named(summary_df, c("table", "n_rows"))
  expect_equal(nrow(summary_df), 3)
  expect_equal(summary_df$n_rows[summary_df$table == "variants"], 3)
  expect_equal(summary_df$n_rows[summary_df$table == "annotations"], 2)
  expect_equal(summary_df$n_rows[summary_df$table == "meta"], 0)

  # Cleanup
  unlink(db_path)
})
