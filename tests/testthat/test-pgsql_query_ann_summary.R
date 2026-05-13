test_that("pgsql_query_ann_summary returns a valid list of summaries", {
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

      # Use recorded coordinates
      res <- pgsql_query_ann_summary(
        con,
        chrom = "Chr03",
        start = 79037800,
        end = 79037890
      )

      # Structural Check
      expect_type(res, "list")
      expect_named(
        res,
        c("annotation_summary", "impact_summary", "variant_type_totals")
      )

      # Content Check for Data Frames
      expect_s3_class(res$annotation_summary, "data.frame")
      expect_s3_class(res$impact_summary, "data.frame")
      expect_s3_class(res$variant_type_totals, "data.frame")

      # Column Check
      expect_true("variant_type" %in% colnames(res$impact_summary))
      expect_true("count" %in% colnames(res$impact_summary))

      DBI::dbDisconnect(con)
    })
  })
})

test_that("pgsql_query_ann_summary handles missing tables correctly", {
  # Logic check using SQLite (No fixtures required)
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # Ensure connection is closed cleanly when test finishes
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Tables don't exist yet, so we expect the dynamic error message
  expect_error(
    pgsql_query_ann_summary(con, chrom = "Chr01", start = 1, end = 100),
    "Table 'annotations' not found"
  )
})