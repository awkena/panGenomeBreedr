test_that("pgsql_summarize_tables returns correct row counts from mock DB", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")

  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      # Connect with recording credentials
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # Execute summary
      summary_df <- pgsql_summarize_tables(con)

      # Validate structure
      expect_s3_class(summary_df, "data.frame")
      expect_equal(colnames(summary_df), c("table", "n_rows"))

      # Check if core tables are present in the summary
      expect_true("variants" %in% summary_df$table)
      expect_true("genotypes" %in% summary_df$table)

      # Validate data types (ensure n_rows is numeric/double to prevent overflow)
      expect_type(summary_df$n_rows, "double")
      expect_true(all(summary_df$n_rows >= 0))

      DBI::dbDisconnect(con)
    })
  })
})