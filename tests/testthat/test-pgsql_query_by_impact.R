test_that("pgsql_query_by_impact filters correctly for every impact level", {
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

      # All levels to test
      target_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

      for (lvl in target_levels) {
        res <- pgsql_query_by_impact(
          con,
          impact_level = lvl,
          chrom = "Chr03",
          start = 79037800,
          end = 79037890
        )

        # 1. Structure check
        expect_s3_class(res, "data.frame")

        # 2. Filtering check
        if (nrow(res) > 0) {
          # Every row in the result MUST match the requested level
          expect_true(
            all(res$impact == lvl),
            info = paste("Failed on level:", lvl)
          )

          # Check that vital columns exist
          expect_true(all(
            c("gene_name", "impact", "variant_id") %in% colnames(res)
          ))
        }
      }

      # 3. Test the default behavior (Full Pangenome Impact)
      res_all <- pgsql_query_by_impact(
        con,
        chrom = "Chr03",
        start = 79037800,
        end = 79037890
      )
      if (nrow(res_all) > 0) {
        expect_true(all(res_all$impact %in% target_levels))
      }

      DBI::dbDisconnect(con)
    })
  })
})