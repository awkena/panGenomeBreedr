test_that("pgsql_query_db retrieves and unpacks genomic data correctly", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")
  skip("Mock SQL needs re-recording")

  # Standard fixture-based test for real data
  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # Target region for tests (recorded fixtures)
      chr <- "Chr03"
      st <- 79037800
      en <- 79037890

      # 1. Test Variants Query
      res_var <- pgsql_query_db(con, "variants", chr, st, en)
      expect_s3_class(res_var, "data.frame")
      expect_true(all(res_var$chrom == chr))

      # 2. Test Annotations Query
      res_ann <- pgsql_query_db(
        con,
        "annotations",
        chr,
        st,
        en,
        gene_name = "Sobic.003G421300"
      )
      expect_true(all(res_ann$gene_name == "Sobic.003G421300"))
      expect_true("impact" %in% colnames(res_ann))

      # 3. Test Genotypes Unpacking (The heavy lifting)
      res_gt <- pgsql_query_db(con, "genotypes", chr, st, en)
      expect_s3_class(res_gt, "data.frame")
      expect_false("calls" %in% colnames(res_gt))
      expect_true(ncol(res_gt) > 1000)

      DBI::dbDisconnect(con)
    })
  })
})

test_that("pgsql_query_db handles empty results gracefully", {
  # We use in-memory SQLite here so we don't need a recorded fixture
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # Create the necessary table structure so the query doesn't crash
  DBI::dbExecute(con, "CREATE TABLE variants (chrom TEXT, pos INTEGER)")

  # Query a chromosome that is logically empty
  res_empty <- pgsql_query_db(con, "variants", "Chr99", 1, 100)

  # Assertions
  expect_s3_class(res_empty, "data.frame")
  expect_equal(nrow(res_empty), 0)

  DBI::dbDisconnect(con)
})
