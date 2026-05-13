test_that("pgsql_variant_impact_summary returns reshaped impact data from mock DB", {
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

      # Execute summary function
      impact_stats <- pgsql_variant_impact_summary(con)

      # 1. Validate structure
      expect_s3_class(impact_stats, "data.frame")
      expect_true("chrom" %in% colnames(impact_stats))

      # 2. Validate wide-format column naming logic
      # Check if any column starts with 'impact_'
      impact_cols <- grep("^impact_", colnames(impact_stats), value = TRUE)
      expect_true(length(impact_cols) > 0)

      # 3. Validate content logic
      if (nrow(impact_stats) > 0) {
        # Ensure no NAs remain (they should be 0)
        expect_false(any(is.na(impact_stats)))
      }

      DBI::dbDisconnect(con)
    })
  })
})