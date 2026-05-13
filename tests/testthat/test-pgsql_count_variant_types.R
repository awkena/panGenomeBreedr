test_that("pgsql_count_variant_types retrieves variant distribution from mock DB", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")

  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # Execute
      res <- pgsql_count_variant_types(con)

      # Structural Check
      expect_s3_class(res, "data.frame")
      expect_equal(colnames(res), c("variant_type", "n"))

      if (nrow(res) > 0) {
        # Ensure counts are numeric (no integer overflow)
        expect_type(res$n, "double")
        # Check for typical types in your pangenome
        expect_true(any(c("SNP", "INDEL") %in% res$variant_type))
      }

      DBI::dbDisconnect(con)
    })
  })
})

test_that("pgsql_count_variant_types handles missing tables and columns", {
  # Logic check using SQLite (No fixtures required)
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # Table doesn't exist
  expect_error(
    pgsql_count_variant_types(con, variants_table = "ghost_table"),
    "does not exist in the database"
  )

  # Test: Column 'variant_type' is missing
  DBI::dbExecute(con, "CREATE TABLE variants (chrom TEXT, pos INTEGER)")
  expect_error(
    pgsql_count_variant_types(con, variants_table = "variants"),
    "does not have a 'variant_type' column"
  )

  DBI::dbDisconnect(con)
})