test_that("pgsql_variant_stats returns correct summary structure from mock DB", {
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

      #  Test full stats 
      stats_all <- pgsql_variant_stats(con, include_annotations = TRUE)

      # Validate structure
      expect_s3_class(stats_all, "data.frame")
      expect_true(all(
        c("chrom", "n_variants", "min_pos", "max_pos", "n_annotated") %in%
          colnames(stats_all)
      ))

      # Validate logic 
      expect_true(nrow(stats_all) > 0)
      if ("Chr03" %in% stats_all$chrom) {
        row_03 <- stats_all[stats_all$chrom == "Chr03", ]
        expect_true(row_03$n_variants > 0)
        expect_true(row_03$max_pos >= row_03$min_pos)
      }

      # Test base stats only 
      stats_base <- pgsql_variant_stats(con, include_annotations = FALSE)

      # Ensure optional column is absent
      expect_false("n_annotated" %in% colnames(stats_base))
      expect_equal(nrow(stats_base), nrow(stats_all))

      DBI::dbDisconnect(con)
    })
  })
})